#ifndef SRC_TOOLS_TIMER
#define SRC_TOOLS_TIMER

#include <chrono>
#include <cassert>

class Timer
{
    bool running = false;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::duration<double> duration = static_cast<std::chrono::duration<double>>(0);

public:
    void start()
    {
        running = true;
        startTime = std::chrono::high_resolution_clock::now();
    }

    void stop()
    {
        const auto endTime = std::chrono::high_resolution_clock::now();
        duration += (endTime - startTime);
        running = false;
    }

    /**
     * Total duration in microseconds.
     */
    [[nodiscard]] auto getTotalDurationUs() const noexcept
    {
        assert(!running);
        return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    }

    [[nodiscard]] auto getTotalDurationS() const noexcept
    {
        assert(!running);
        return duration.count();
    }
};

class ScopedTimer
{
    Timer& timer;
public:
    explicit ScopedTimer(Timer& timer)
        : timer(timer)
    {
        timer.start();
    }

    ~ScopedTimer()
    {
        timer.stop();
    }
};

inline std::ostream& operator<<(std::ostream& stream, const Timer& timer) noexcept
{
    return stream << timer.getTotalDurationS() << "s (" << timer.getTotalDurationUs() << "us)";
}


#endif //SRC_TOOLS_TIMER
