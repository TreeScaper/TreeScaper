#pragma once

#include <string>
#include <curl/curl.h>

using namespace std;

extern const string bootstrap_results_filename;

class CRAHandle {
public:
	CRAHandle(string appkey, string username, string password) :
		appkey(appkey),
		username(username),
		password(password)
		{};
	void submit_raxml(string filename);
	void retrieve_raxml(string jobname);
private:
	CURL *create_cra_curl_handle();
	string appkey;
	string username;
	string password;
	string userdata;
};
