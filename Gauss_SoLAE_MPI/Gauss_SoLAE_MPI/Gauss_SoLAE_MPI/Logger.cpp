#define _CRT_SECURE_NO_WARNINGS
#include<ctime>
#include "Logger.h"

using namespace std;

string Logger::path;

Logger::Logger() {
	init(
		path.empty()
		? "logfile.txt"
		: string(path).append("logfile.txt").data()
	);
}

Logger::Logger(const char *filename) {
	init(filename);
}

Logger::Logger(string path) {
	Logger::path = path;
	init(
		path.empty()
		? "logfile.txt"
		: string(path).append("logfile.txt").data()
	);
}

void Logger::init(const char* filename) {
	try {
		if (!_fhnd.is_open())
		{
			_fhnd.open(filename, ofstream::out | ofstream::app);
			strcpy(_dt_fmt, "[%Y-%m-%d %H:%M:%S]");
		}
	}
	catch (const std::exception& e)
	{
		cout <<"Exception!!! " << endl;
	}
}

Logger::~Logger() {
	_fhnd.close();
}

void Logger::setDateFormat(char *fmt) {
	strcpy(_dt_fmt, fmt);
}

//GetTimeStamp obtains the current timestamp and formats to "[2014-09-28 22:12:59]" form
string Logger::getTimeStamp() {
	char buff[100];
	time_t t = time(nullptr);
	tm tm;
	localtime_s(&tm, &t);
	strftime(buff, sizeof(buff), _dt_fmt, &tm);
	return string(buff);
}

string Logger::getLoglevelString(loglevel level) {
	string lvlstr{};
	switch (level) {
	case LOGWARN: lvlstr = "WARN: "; break;
	case LOGERROR: lvlstr = "ERROR: "; break;
	case LOGINFO: lvlstr = "INFO: "; break;
	}
	return lvlstr;
}

void Logger::log(loglevel level, string message) {
	_fhnd << getTimeStamp() << " " << getLoglevelString(level) << message << endl;
}

void Logger::logI(string message) {
	log(LOGINFO, message);
}

void Logger::logE(string message) {
	log(LOGERROR, message);
}

void Logger::logW(string message) {
	log(LOGWARN, message);
}