#ifndef TIMER_H_
#define TIMER_H_

#ifndef __PGI
#include <atomic>
#endif

#include <chrono>

namespace rdmini {
namespace timer {

inline void fence() {
#ifndef __PGI
    std::atomic_thread_fence(std::memory_order_seq_cst);
#else
    __sync_synchronize();
#endif
}

struct hr_timer {
    typedef double value_type;

    hr_timer() { reset(); }
    
    /** Reset accumulator, stop counting */
    void reset() {
        fence();
        paused=true;
        zero_accumulator();
        fence();
    }

    /** Resume counting from paused state */
    void resume() {
        fence();
        paused=false;
        read_counter(c0);
        fence();
    }

    /** Pause counting */
    void stop() {
        fence();
        read_counter(c1);
        a+=c1-c0;
        paused=true;
        fence();
    }

    /** Return accumulated counter difference
     *
     * Only valid when paused.
     */
    value_type time() {
        fence();
        return convert(a);
    }

    void start() {
        reset();
        resume();
    }


private:
    typedef std::chrono::high_resolution_clock clock;
    typedef clock::time_point count_type;
    typedef clock::duration accum_type;

    count_type c0,c1;
    accum_type a;
    bool paused;

    value_type convert(accum_type a) {
        return std::chrono::duration_cast<std::chrono::duration<value_type>>(a).count();
    }

    void zero_accumulator() {
        a=accum_type::zero();
    }

    void read_counter(count_type &c) {
        c=clock::now();
    }
};

template <typename Timer>
struct timer_guard {
    Timer *Tp;

    timer_guard(Timer &T_): Tp(&T_) {
        fence();
        Tp->start();
        fence();
    }

    timer_guard(timer_guard &&g): Tp(g.Tp) { g.Tp=nullptr; }
    timer_guard &operator=(timer_guard &&g) { if (this!=&g) g.Tp=nullptr; return *this; }

    timer_guard() =delete;
    timer_guard(const timer_guard &g) =delete;
    timer_guard &operator=(const timer_guard &g) =delete;

    ~timer_guard() {
        if (Tp) Tp->stop();
    }
};

template <typename Timer>
inline timer_guard<Timer> guard(Timer &T) { return timer_guard<Timer>(T); }

}} // namespace rdmini::timer


#endif // ndef TIMER_H_

