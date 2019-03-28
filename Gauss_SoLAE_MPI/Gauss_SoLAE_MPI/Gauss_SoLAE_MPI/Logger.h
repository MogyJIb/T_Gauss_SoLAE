#pragma once
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<mutex>

using namespace std;

class Logger
{

	static string path;

public:
	enum loglevel { LOGWARN, LOGERROR, LOGINFO };
	ofstream _fhnd;

	void setDateFormat(char*);
	void logI(string);
	void logE(string);
	void logW(string);

	Logger(string path);
	Logger();
	Logger(const char*);
	~Logger();

private:
	char _dt_fmt[512];
	const std::vector<string>_loglevel_text;
	string getTimeStamp();
	string getLoglevelString(loglevel);
	void log(loglevel, string);
	void init(const char*);
};

