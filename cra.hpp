#pragma once

#include <string>
#include <chrono>
#include <curl/curl.h>
#include <pugixml.hpp>

using namespace std;

// TreeScaper's classification of job status.
enum JobStatus {
	UNSUBMITTED,
	SUBMITTED,
	SUCCESSFUL,
	FAILED
};

enum CRALogLevel {
	NONE,
	DEBUG
};

extern enum CRALogLevel cra_log_level;

// Stores information about a CRA job.
class CRAJob {
	public:
		CRAJob(string inputfile, enum JobStatus status = UNSUBMITTED, string handle = "") :
			inputfile(inputfile),
			status(status),
			handle(handle) {};

		// Input file for CRA job.
		// Most jobs have a single input file.
		string inputfile;

		// Status of job as tracked in TreeScaper.
		enum JobStatus status;

		// Handle for constructing url to find status.
		string handle;
};

// Manages interaction with the CRA.
class CRAHandle {
public:
	// Initializes authentication and CRA application ID.
	CRAHandle() :
		active_jobs(0),
		min_poll_interval_seconds(60)
		{
			// Get CRA username supplied as environment variable.
			const char *username_env = getenv("TS_CRA_USERNAME");
			if (username_env == NULL) {
				throw invalid_argument("Must supply TS_CRA_USERNAME environment variable.");
			}
			username = string(username_env);

			// Get CRA password supplied as environment variable.
			const char *password_env = getenv("TS_CRA_PASSWORD");
			if (password_env == NULL) {
				throw invalid_argument("Must supply TS_CRA_PASSWORD environment variable.");
			}
			password = string(password_env);
		};
	bool submit_jobs(string filelist, string paramfile);
private:
	bool submit_job(CRAJob& job);
	bool parse_single_status(pugi::xml_node status_node);
	bool parse_status_list();
	bool retrieve_url(string url);
	bool change_job_status(CRAJob& job, enum JobStatus status);
	bool write_job_status();
	CURL *create_cra_curl_handle();

	// Username provided by user.
	string username;

	// Password provided by user.
	string password;

	// Buffer to store response data from CURL request.
	string userdata;

	// Buffer to store response headers from CURL request.
	string headerdata;

	// Number of CRA jobs currently active.
	// This must stay below some threshold.
	int active_jobs;

	// List of jobs to run.
	vector<CRAJob> jobs;

	// Parameters for tool
	map<string, string> cra_params;

	// Minimum number of seconds between polls to check status of one job
	int min_poll_interval_seconds;

	// Last time status of all active jobs in the jobs list was checked.
	// We rate limit this check.
	chrono::time_point<std::chrono::system_clock> last_poll;

	// Logs last HTTP response data and header.
	void log_last_response();
};
