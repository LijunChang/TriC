#include "Timer.h"

Timer::Timer() {
	m_start = timestamp();
}

void Timer::restart() {
	m_start = timestamp();
}

long long Timer::elapsed() {
	return timestamp() - m_start;
}

long long Timer::timestamp() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
}
