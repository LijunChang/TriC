/******************************************************************************
 * timer.h 
 * record the elapsed time in microseconds (10^{-6} second)
 *****************************************************************************/

#ifndef _TIMER_LJ_
#define _TIMER_LJ_

#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
	Timer() ;
	void restart() ;
	long long elapsed() ;

private:
	long long m_start;

	// Returns a timestamp ('now') in microseconds
	long long timestamp() ;
};

#endif
