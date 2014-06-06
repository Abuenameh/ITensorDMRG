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

typedef function<void(ostream&, istream&)> callback;
//typedef stream<file_descriptor_sink> postream;
//typedef stream<file_descriptor_source> pistream;

class ProcessPool {
public:
    ProcessPool(size_t, string, callback);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> future<typename result_of<F(ostream&, istream&, Args...)>::type>;
    void interrupt();
    ~ProcessPool();
private:
    vector<child> children;
    vector<ostream> oss;
    vector<istream> iss;
    //vector<file_descriptor_sink> ofds;
    //vector<file_descriptor_source> ifds;
    //vector<outfd> ofds;
    //vector<infd> ifds;
    vector<int> ofds;
    vector<int> ifds;
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
    for(size_t i = 0;i<processes;++i) {
        
        boost::process::pipe pin = create_pipe();
        file_descriptor_sink insink(pin.sink, never_close_handle);
        file_descriptor_source insource(pin.source, never_close_handle);
        
        boost::process::pipe pout = create_pipe();
        file_descriptor_sink outsink(pout.sink, never_close_handle);
        file_descriptor_source outsource(pout.source, never_close_handle);
        
        
        cout << "Sink: " << pin.sink << endl << "Source: " << pout.source << endl;
        
        stream<file_descriptor_sink> os(insink);
        //oss.push_back(os);
        //oss.emplace_back(insink);
        
        stream<file_descriptor_source> is(outsource);
        //iss.push_back(is);
        //iss.emplace_back(outsource);
        
        //ofds.push_back(insink);
        //ifds.push_back(outsource);
        
        ofds.push_back(pin.sink);
        ifds.push_back(pout.source);
        
        children.push_back(execute(set_args(args), inherit_env(), bind_stdin(insource), bind_stdout(outsink)));
        
        init(os, is);
        
        workers.emplace_back(
            [this,i]
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
                    stream<outfd> os(outfd(this->ofds[i], never_close_handle));
                    stream<infd> is(infd(this->ifds[i], never_close_handle));
                    task(os, is);
                    //task(this->ofds[i], this->ifds[])
                }
            }
        );
    }
}

// add new work item to the pool
template<class F, class... Args>
auto ProcessPool::enqueue(F&& f, Args&&... args) 
    -> future<typename result_of<F(ostream&, istream&, Args...)>::type>
{
    typedef typename result_of<F(ostream&, istream&, Args...)>::type return_type;
    
    // don't allow enqueueing after stopping the pool
    if(stop)
        throw runtime_error("enqueue on stopped ProcessPool");

    auto task = std::make_shared< std::packaged_task<return_type(ostream&, istream&)> >(
            bind(forward<F>(f), placeholders::_1, placeholders::_2, forward<Args>(args)...)
        );
        
    future<return_type> res = task->get_future();
    {
        unique_lock<mutex> lock(queue_mutex);
        tasks.push([task](ostream& os, istream& is){ (*task)(os, is); });
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