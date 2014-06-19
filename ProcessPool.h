#ifndef PROCESS_POOL_H
#define PROCESS_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

#include <boost/process.hpp>
#include <boost/iostreams/stream.hpp>

#undef SP
#undef TCP_NODELAY
#include <nnxx/socket>

using namespace std;
using namespace std::placeholders;
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

//using namespace nnxx;

typedef file_descriptor_sink outfd;
typedef file_descriptor_source infd;
//typedef function<void(outfd&, infd&)> callback;

//typedef function<void(ostream&, istream&, istream&, bool&)> callback;
//typedef function<void(nnxx::socket&, nnxx::socket&, bool&)> callback;
typedef function<void(message_queue&, bool&)> callback;
//typedef function<void(socket_t&, socket_t&, bool&)> callback;
//typedef stream<file_descriptor_sink> postream;
//typedef stream<file_descriptor_source> pistream;

class ProcessPool
{
public:
    ProcessPool(size_t, string, /*vector<socket_t>&, vector<context_t>&, vector<socket_t>&,*/ callback);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
    -> future<typename result_of<F(message_queue&, bool&, Args...)>::type>;
//    -> future<typename result_of<F(nnxx::socket&, nnxx::socket&, bool&, Args...)>::type>;
    void create_process(string prog, callback init);
    void interrupt();
    ~ProcessPool();

    static int queue_idx;
private:
    vector<child> children;
    // need to keep track of threads so we can join them
    vector<thread> workers;
    // the task queue
    queue<callback> tasks;
    
    vector<context_t> outcontexts;
    vector<socket_t> outsockets;
    vector<context_t> incontexts;
    vector<socket_t> insockets;
    vector<connection_monitor> monitors;
    vector<thread> monitorthreads;
    
    vector<nnxx::socket> osockets;
    vector<nnxx::socket> isockets;
    
    vector<unique_ptr<message_queue> > queues;
    vector<string> queue_names;

    // synchronization
    mutex queue_mutex;
    condition_variable condition;
    bool stop;
    
    static int port;
};

int ProcessPool::port = 5557;

int ProcessPool::queue_idx = 0;

void ProcessPool::create_process(string prog, callback init)
{
    string outport = to_string(port++);
    string inport = to_string(port++);
    
    osockets.emplace_back(nnxx::SP, nnxx::PUSH);
    isockets.emplace_back(nnxx::SP, nnxx::PULL);
    
    nnxx::socket& os = osockets.back();
    nnxx::socket& is = isockets.back();
    
    os.connect("tcp://127.0.0.1:" + outport);
    is.bind("tcp://127.0.0.1:" + inport);
    
    string queue_name = "dmrg." + to_string(queue_idx++);
    queue_names.push_back(queue_name);
//    queues.emplace_back(create_only, queue_name.c_str(), 100, MAX_MSG_SIZE);
    unique_ptr<message_queue> queue = unique_ptr<message_queue>(new message_queue(create_only, queue_name.c_str(), 1000, MAX_MSG_SIZE));
    message_queue& mq = *queue;
    queues.push_back(move(queue));
//    message_queue& mq = queues.back();

    /*outcontexts.emplace_back();
    outsockets.emplace_back(outcontexts.back(), ZMQ_PUSH);
        
    incontexts.emplace_back();
    insockets.emplace_back(incontexts.back(), ZMQ_PULL);
        
    context_t& ctx = incontexts.back();
    string monitoraddr = "inproc://monitor." + inport;
    //monitors.emplace_back(incontexts.back(), insockets.back());
    monitors.emplace_back([&] () {
        cout << "---Disconnected---" << endl << flush;
        ctx.close();
    });
        monitorthreads.emplace_back([&] () {
            monitors.back().monitor(insockets.back(), monitoraddr.c_str(), ZMQ_EVENT_DISCONNECTED);
        });
    
        socket_t& os = outsockets.back();
        os.connect(("tcp://127.0.0.1:" + outport).c_str());
        
        socket_t& is = insockets.back();
        is.bind(("tcp://127.0.0.1:" + inport).c_str());*/
    
    /*int count = 0;
    cout << "Creating process" << endl;*/
#ifdef FST
    children.push_back(execute(set_args(vector<string> {prog, queue_name, outport, inport}), inherit_env()/*, bind_stdin(file_descriptor_source(outpipe.source, never_close_handle)), bind_stdout(file_descriptor_sink(inpipe.sink, never_close_handle))*/));
#else
    children.push_back(execute(set_args(vector<string> {prog}), inherit_env(), bind_stdin(file_descriptor_source(outpipe.source, never_close_handle)), bind_stdout(file_descriptor_sink(inpipe.sink, never_close_handle)), bind_fd(42, file_descriptor_sink(abortpipe.sink, never_close_handle))));
#endif
    child& c = children.back();

    /*cout << "Here " << count++ << endl;
    bool abort;
    init(os, is, abortis, abort);*/
    
    bool abort;
    init(mq, abort);

    workers.emplace_back(
    [this,prog,init,&mq/*&os,&is*//*outpipe,inpipe,abortpipe*/] {
        for(;;) {
            unique_lock<mutex> lock(this->queue_mutex);
            while(!this->stop && this->tasks.empty())
                this->condition.wait(lock);
            if(this->stop || this->tasks.empty())
                return;
            callback task(this->tasks.front());
            cout << "Got task" << endl << flush;
            this->tasks.pop();
            lock.unlock();

            /*stream<file_descriptor_sink> os(file_descriptor_sink(outpipe.sink, never_close_handle));
            stream<file_descriptor_source> is(file_descriptor_source(inpipe.source, never_close_handle));
            stream<file_descriptor_source> abortis(file_descriptor_source(abortpipe.source, never_close_handle));

            bool abort;
            task(os, is, abortis, abort);
            if(abort) {
                lock.lock();
                this->create_process(prog, init);
                this->tasks.push(task);
                lock.unlock();
                break;
            }*/
            
//            try {
//                cout << "Trying" << endl;
            bool abort = false;
            cout << "Running task " << (int)(&mq) << endl << flush;
            task(mq, abort);
            cout << "Aborted " << (int)(&mq) << ": " << abort << endl << flush;
            if(abort && !this->stop) {
                lock.lock();
                cout << "Recreating process" << endl;
                this->create_process(prog, init);
                this->tasks.push(task);
                lock.unlock();
                break;
            }
//            cout << "Tried" << endl;
//            } catch (...) {cout << "Exception" << endl;}
        }
    }
    );
}


// the constructor just launches some amount of workers
inline ProcessPool::ProcessPool(size_t processes, string prog, /*vector<socket_t>& outsockets, vector<context_t>& incontexts, vector<socket_t>& insockets,*/ callback init)
    :   /*outsockets(outsockets), incontexts(incontexts), insockets(insockets),*/ stop(false)
{
    for(size_t i = 0; i<processes; ++i) {
        cout << "Creating initial process" << endl;
        create_process(prog, init);
    }
}

// add new work item to the pool
template<class F, class... Args>
auto ProcessPool::enqueue(F&& f, Args&&... args)
-> future<typename result_of<F(message_queue&, /*nnxx::socket&, nnxx::socket&,*/ bool&, Args...)>::type> {
    typedef typename result_of<F(message_queue&, /*nnxx::socket&, nnxx::socket&,*/ bool&, Args...)>::type return_type;

    cout << "Enqueueing" << endl << flush;
    // don't allow enqueueing after stopping the pool
    if(stop)
        throw runtime_error("enqueue on stopped ProcessPool");

    auto task = std::make_shared< std::packaged_task<return_type(message_queue&, /*nnxx::socket&, nnxx::socket&,*/ bool&)> >(
        bind(forward<F>(f), placeholders::_1, placeholders::_2, forward<Args>(args)...)
//        bind(forward<F>(f), placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4, forward<Args>(args)...)
    );

    future<return_type> res = task->get_future();
    {
        unique_lock<mutex> lock(queue_mutex);
        tasks.push([task,&res](message_queue& mq, /*nnxx::socket& os, nnxx::socket& is,*/ bool& abort) {
//            try {
            (*task)(mq, abort);
//            } catch(...) { cout<<"Exception"<<endl;}
            if(abort) {
                task->reset();
            }
        });
    }
    condition.notify_one();
    return res;
}

void ProcessPool::interrupt()
{
    {
        unique_lock<mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(size_t i = 0; i < children.size(); ++i) {
        try {
            terminate(children[i]);
            wait_for_exit(children[i]);
        } catch(...) {}
        message_queue::remove(queue_names[i].c_str());
    }
    //for(size_t i = 0; i<workers.size(); ++i)
        //workers[i].join();
}

// the destructor joins all threads
inline ProcessPool::~ProcessPool()
{
    {
        unique_lock<mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(size_t i = 0; i < children.size(); ++i) {
        try {
            wait_for_exit(children[i]);
        } catch(...) {}
        message_queue::remove(queue_names[i].c_str());
    }
    //for(size_t i = 0; i<workers.size(); ++i)
        //workers[i].join();
}

#endif
