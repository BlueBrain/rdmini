#include <iostream>
#include <thread>

#include "rdmini/timer.h"

namespace timer=rdmini::timer;

void sleep(double dt) {
    std::this_thread::sleep_for(std::chrono::duration<double>(dt));
}

int main() {
    timer::hr_timer T;

    std::cout << "T.time()=" << T.time() << "\n";
    T.start();
    sleep(2.3);
    T.stop();
    std::cout << "After sleep(2.3): T.time()=" << T.time() << "\n";
    T.resume();
    sleep(3.5);
    T.stop();
    std::cout << "After sleep(3.5): T.time()=" << T.time() << "\n";
    T.reset();
    std::cout << "Reset timer.\n";
    {
        auto _=timer::guard(T);
        sleep(4.8);
    }
    std::cout << "After sleep(4.8) (using guard): T.time()=" << T.time() << "\n";
}

