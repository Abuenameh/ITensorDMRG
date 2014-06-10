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

using namespace std;
using namespace std::placeholders;
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

typedef file_descriptor_sink outfd;
typedef file_descriptor_source infd;
//typedef function<void(outfd&, infd&)> callback;

typedef function<void(ostream&, istream&, istream&, bool&)> callback;
//typedef stream<file_descriptor_sink> postream;
//typedef stream<file_descriptor_source> pistream;

class ProcessPool
{
public:
    ProcessPool(size_t, string, callback);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
    -> future<typename result_of<F(ostream&, istream&, istream&, bool&, Args...)>::type>;
    void create_process(string prog, callback init);
    void interrupt();
    ~ProcessPool();
private:
    vector<child> children;
    // need to keep track of threads so we can join them
    vector<thread> workers;
    // the task queue
    queue<callback> tasks;

    // synchronization
    mutex queue_mutex;
    condition_variable condition;
    bool stop;
};

void ProcessPool::create_process(string prog, callback init)
{
    boost::process::pipe outpipe = create_pipe();
    stream<file_descriptor_sink> os(file_descriptor_sink(outpipe.sink, never_close_handle));

    boost::process::pipe inpipe = create_pipe();
    stream<file_descriptor_source> is(file_descriptor_source(inpipe.source, never_close_handle));

    boost::process::pipe abortpipe = create_pipe();
    stream<file_descriptor_source> abortis(file_descriptor_source(abortpipe.source, never_close_handle));

    children.push_back(execute(set_args(vector<string> {prog}), inherit_env(), bind_stdin(file_descriptor_source(outpipe.source, never_close_handle)), bind_stdout(file_descriptor_sink(inpipe.sink, never_close_handle)), bind_fd(42, file_descriptor_sink(abortpipe.sink, never_close_handle))));

    bool abort;
    init(os, is, abortis, abort);

    workers.emplace_back(
    [this,prog,init,outpipe,inpipe,abortpipe] {
        for(;;) {
            unique_lock<mutex> lock(this->queue_mutex);
            while(!this->stop && this->tasks.empty())
                this->condition.wait(lock);
            if(this->stop || this->tasks.empty())
                return;
            callback task(this->tasks.front());
            this->tasks.pop();
            lock.unlock();

            stream<file_descriptor_sink> os(file_descriptor_sink(outpipe.sink, never_close_handle));
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
            }
        }
    }
    );
}


// the constructor just launches some amount of workers
inline ProcessPool::ProcessPool(size_t processes, string prog, callback init)
    :   stop(false)
{
    for(size_t i = 0; i<processes; ++i) {
        create_process(prog, init);
    }
}

// add new work item to the pool
template<class F, class... Args>
auto ProcessPool::enqueue(F&& f, Args&&... args)
-> future<typename result_of<F(ostream&, istream&, istream&, bool&, Args...)>::type> {
    typedef typename result_of<F(ostream&, istream&, istream&, bool&, Args...)>::type return_type;

    // don't allow enqueueing after stopping the pool
    if(stop)
        throw runtime_error("enqueue on stopped ProcessPool");

    auto task = std::make_shared< std::packaged_task<return_type(ostream&, istream&, istream&, bool&)> >(
        bind(forward<F>(f), placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4, forward<Args>(args)...)
    );

    future<return_type> res = task->get_future();
    {
        unique_lock<mutex> lock(queue_mutex);
        tasks.push([task,&res](ostream& os, istream& is, istream& abortis, bool& abort) {
            (*task)(os, is, abortis, abort);
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
    }
    //for(size_t i = 0; i<workers.size(); ++i)
        //workers[i].join();
}

#endif
