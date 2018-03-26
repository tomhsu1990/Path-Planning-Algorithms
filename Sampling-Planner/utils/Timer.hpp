#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <time.h>

class Timer {
public:
	Timer(){
		this->start();
	}
	virtual ~Timer(){};

	void start() {
		m_start_t = clock();
		gettimeofday(&m_start_time, NULL);
	}

	void restart() {
		this->start();
	}

	double getElapsedMilliseconds(){
		timeval now;
		gettimeofday(&now, NULL);
		double milliseconds = ((now.tv_sec - m_start_time.tv_sec)*1000000.0 + (now.tv_usec - m_start_time.tv_usec))/1000.0;
		return milliseconds;
	}

	double getElapsedCPUMilliseconds(){
		clock_t now = clock();
		double milliseconds = 1000.0*(now - m_start_t) / CLOCKS_PER_SEC;
		return milliseconds;
	}
private:
	clock_t m_start_t;
	timeval m_start_time;
};

#endif // TIMER_H_
