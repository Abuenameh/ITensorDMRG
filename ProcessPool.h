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

#include <zmq.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/process.hpp>

using namespace std;
using namespace zmq;
using namespace boost::iostreams;
using namespace boost::process;
using namespace boost::process::initializers;

typedef function<void(socket_t&, socket_t&, ostream&)> callback;

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
    vector<socket_t> socketsout;
    vector<socket_t> socketsin;
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
    context_t context;
    string addr;
    //socket_t socketin;
    //socket_t socketout;
    int port = 5555;
    vector< string > args;
    args.push_back(prog);
    args.push_back("");
    for(size_t i = 0;i<processes;++i) {
        args[1] = to_string(port + 2*i);
        
        /*socketin = socket_t(context, ZMQ_PULL);
        addr = "tcp://localhost:" + to_string(port + 2*i);
        socketin.connect(addr.c_str());
        socketsin.emplace_back(socketin);*/
        
        /*socketout = socket_t(context, ZMQ_PUSH);
        addr = "tcp://localhost:" + to_string(port + 2*i + 1);
        socketin.connect(addr.c_str());
        socketsout.emplace_back(socketout);*/
        
        socketsout.emplace_back(context, ZMQ_PUSH);
        addr = "tcp://localhost:" + to_string(port + 2*i);
        //socketsout.back().connect(addr.c_str());

        socketsin.emplace_back(context, ZMQ_PULL);
        addr = "tcp://localhost:" + to_string(port + 2*i + 1);
        //socketsin.back().bind(addr.c_str());
        
        boost::process::pipe p = create_pipe();
        file_descriptor_sink sink(p.sink, never_close_handle);
        
        file_descriptor_source source(p.source, never_close_handle);
        
        stream<file_descriptor_sink> os(sink);
        
        children.push_back(execute(set_args(args), inherit_env(), bind_stdin(source)));
        
        init(socketsout.back(), socketsin.back(), os);
        
        workers.emplace_back(
            [this,i,&os]
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
                    task(this->socketsout[i], this->socketsin[i], os);
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