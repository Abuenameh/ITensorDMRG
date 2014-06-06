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

//#include <zmq.hpp>

#include <boost/process.hpp>
#include <boost/iostreams/stream.hpp>

using namespace std;
using namespace zmq;
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

typedef function<void(ostream&, istream&)> callback;
typedef stream<file_descriptor_sink> postream;
typedef stream<file_descriptor_source> pistream;

class ProcessPool {
public:
    ProcessPool(size_t, string, callback);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> future<typename result_of<F(Args...)>::type>;
    void interrupt();
    ~ProcessPool();
private:
    vector<child> children;
    vector<ostream> oss;
    vector<istream> iss;
    // need to keep track of threads so we can join them
    vector<thread> workers;
    // the task queue
    queue<callback> tasks;
    
    // synchronization
    mutex queue_mutex;
    condition_variable condition;
    bool stop;
};
 
// the constructor just launches some amount of workers
inline ProcessPool::ProcessPool(size_t processes, string prog, callback init)
    :   stop(false)
{
    vector< string > args;
    args.push_back(prog);
//    args.push_back("");
    for(size_t i = 0;i<processes;++i) {
//        args[1] = to_string(port + 2*i);
        
        boost::process::pipe pin = create_pipe();
        file_descriptor_sink insink(pin.sink, never_close_handle);
        file_descriptor_source insource(pin.source, never_close_handle);
        
        boost::process::pipe pout = create_pipe();
        file_descriptor_sink outsink(pout.sink, never_close_handle);
        file_descriptor_source outsource(pout.source, never_close_handle);
        
        stream<file_descriptor_sink> os(insink);
        //oss.push_back(os);
        //oss.emplace_back(insink);
        
        stream<file_descriptor_source> is(outsource);
        //iss.push_back(is);
        //iss.emplace_back(outsource);
        
        children.push_back(execute(set_args(args), inherit_env(), bind_stdin(insource), bind_stdout(outsink)));
        
        //init(oss.back(), iss.back());
        init(os, is);
        
        wait_for_exit(children.back());
        
        workers.emplace_back(
            [this,i,&os,&is]
            {
                for(;;)
                {
                    unique_lock<mutex> lock(this->queue_mutex);
                    while(!this->stop && this->tasks.empty())
                        this->condition.wait(lock);
                    if(this->stop || this->tasks.empty())
                        return;
                    callback task(this->tasks.front());
                    this->tasks.pop();
                    lock.unlock();
                    //task(this->oss[i], this->iss[i]);
                    task(os, is);
                }
            }
        );
    }
}

// add new work item to the pool
template<class F, class... Args>
auto ProcessPool::enqueue(F&& f, Args&&... args) 
    -> future<typename result_of<F(Args...)>::type>
{
    typedef typename result_of<F(Args...)>::type return_type;
    
    // don't allow enqueueing after stopping the pool
    if(stop)
        throw runtime_error("enqueue on stopped ProcessPool");

    auto task = make_shared< packaged_task<return_type()> >(
            bind(forward<F>(f), forward<Args>(args)...)
        );
        
    future<return_type> res = task->get_future();
    {
        unique_lock<mutex> lock(queue_mutex);
        tasks.push([task](){ (*task)(); });
    }
    condition.notify_one();
    return res;
}

void ProcessPool::interrupt()
{
    unique_lock<mutex> lock(queue_mutex);
    stop = true;
}

// the destructor joins all threads
inline ProcessPool::~ProcessPool()
{
    {
        unique_lock<mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(size_t i = 0; i < children.size(); ++i)
        wait_for_exit(children[i]);
    for(size_t i = 0;i<workers.size();++i)
        workers[i].join();
}

#endif