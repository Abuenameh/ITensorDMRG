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
#include <chrono>

#include <boost/process.hpp>
#include <boost/iostreams/stream.hpp>

using namespace std;
using namespace std::chrono;
using namespace std::placeholders;
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

typedef function<void(message_queue&, message_queue&, bool&)> callback;
typedef function<void(message_queue&, message_queue&, message_queue&, bool&)> init_callback;


class ProcessPool
{
public:
    ProcessPool(size_t, string, init_callback);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
    -> future<typename result_of<F(message_queue&, message_queue&, bool&, Args...)>::type>;
    void create_process(string prog, init_callback init);
    void interrupt();
    ~ProcessPool();

    static system_clock::rep queue_idx;
//    static system_clock::rep setup_queue_idx;
    static string setup_queue_name;
    
    static int L;
    
private:
    vector<child> children;
    // need to keep track of threads so we can join them
    vector<thread> workers;
    // the task queue
    queue<callback> tasks;
    
    vector<unique_ptr<message_queue> > cqueues;
    vector<unique_ptr<message_queue> > iqueues;
    vector<unique_ptr<message_queue> > oqueues;
    vector<string> cqueue_names;
    vector<string> iqueue_names;
    vector<string> oqueue_names;

    // synchronization
    mutex queue_mutex;
    condition_variable condition;
    bool stop;
    
//    static int port;
};

//int ProcessPool::port = 5557;

system_clock::rep ProcessPool::queue_idx = 0;

string ProcessPool::setup_queue_name = "";

int ProcessPool::L = 0;

void ProcessPool::create_process(string prog, init_callback init)
{
    string cqueue_name = "dmrg." + to_string(queue_idx++);
    message_queue cq(create_only, cqueue_name.c_str(), L, CDS_MAX_MSG_SIZE);
//    cqueue_names.push_back(cqueue_name);
//    unique_ptr<message_queue> cqueue = unique_ptr<message_queue>(new message_queue(create_only, cqueue_name.c_str(), L, CDS_MAX_MSG_SIZE));
//    message_queue& cq = *cqueue;
//    cqueues.push_back(move(cqueue));

    string oqueue_name = "dmrg." + to_string(queue_idx++);
    oqueue_names.push_back(oqueue_name);
    unique_ptr<message_queue> oqueue = unique_ptr<message_queue>(new message_queue(create_only, oqueue_name.c_str(), NUM_MSG, MAX_MSG_SIZE));
    message_queue& oq = *oqueue;
    oqueues.push_back(move(oqueue));

    string iqueue_name = "dmrg." + to_string(queue_idx++);
    iqueue_names.push_back(iqueue_name);
    unique_ptr<message_queue> iqueue = unique_ptr<message_queue>(new message_queue(create_only, iqueue_name.c_str(), NUM_MSG, MAX_MSG_SIZE));
    message_queue& iq = *iqueue;
    iqueues.push_back(move(iqueue));

//#ifdef FST
    children.push_back(execute(set_args(vector<string> {prog, cqueue_name, oqueue_name, iqueue_name}), inherit_env()));
//#else
//    children.push_back(execute(set_args(vector<string> {prog}), inherit_env(), bind_stdin(file_descriptor_source(outpipe.source, never_close_handle)), bind_stdout(file_descriptor_sink(inpipe.sink, never_close_handle)), bind_fd(42, file_descriptor_sink(abortpipe.sink, never_close_handle))));
//#endif
    child& c = children.back();

    bool abort = false;
    init(cq, oq, iq, abort);
    
    message_queue::remove(cqueue_name.c_str());

    workers.emplace_back(
    [this, prog, init, &oq, &iq] {
        for(;;) {
            unique_lock<mutex> lock(this->queue_mutex);
            while(!this->stop && this->tasks.empty())
                this->condition.wait(lock);
            if(this->stop || this->tasks.empty())
                return;
            callback task(this->tasks.front());
            this->tasks.pop();
            lock.unlock();

            bool abort = false;
            task(oq, iq, abort);
            if(abort && !this->stop) {
                lock.lock();
                this->create_process(prog, init);
                this->tasks.push(task);
                lock.unlock();
                break;
            }
        }
    }
    );
}


// the constructor just launches some amount of workers
inline ProcessPool::ProcessPool(size_t processes, string prog, init_callback init)
    :   stop(false)
{
    for(size_t i = 0; i<processes; ++i) {
        create_process(prog, init);
    }
}

// add new work item to the pool
template<class F, class... Args>
auto ProcessPool::enqueue(F&& f, Args&&... args)
-> future<typename result_of<F(message_queue&, message_queue&, bool&, Args...)>::type> {
    typedef typename result_of<F(message_queue&, message_queue&, bool&, Args...)>::type return_type;

    // don't allow enqueueing after stopping the pool
    if(stop)
        throw runtime_error("enqueue on stopped ProcessPool");

    auto task = std::make_shared< std::packaged_task<return_type(message_queue&, message_queue&, bool&)> >(
        bind(forward<F>(f), placeholders::_1, placeholders::_2, placeholders::_3, forward<Args>(args)...)
    );

    future<return_type> res = task->get_future();
    {
        unique_lock<mutex> lock(queue_mutex);
        tasks.push([task,&res](message_queue& oq, message_queue& iq, bool& abort) {
            (*task)(oq, iq, abort);
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
        message_queue::remove(cqueue_names[i].c_str());
        message_queue::remove(oqueue_names[i].c_str());
        message_queue::remove(iqueue_names[i].c_str());
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
        message_queue::remove(cqueue_names[i].c_str());
        message_queue::remove(oqueue_names[i].c_str());
        message_queue::remove(iqueue_names[i].c_str());
    }
    //for(size_t i = 0; i<workers.size(); ++i)
        //workers[i].join();
}

#endif
