#pragma once

#include <string>
#include <chrono>
#include <curl/curl.h>

using namespace std;

// Filename we expect to find RAxML bootstrap results in.
extern const string bootstrap_results_filename;

// TreeScaper's classification of job status.
enum JobStatus {
	UNSUBMITTED,
	SUBMITTED,
	SUCCESSFUL,
	FAILED
};

// Stores information about a CRA job.
class CRAJob {
	public:
		CRAJob(string inputfile) : inputfile(inputfile) {
			status = UNSUBMITTED;
		};

		// Input file for CRA job.
		// Most jobs have a single input file.
		string inputfile;

		// Status of job as tracked in TreeScaper.
		enum JobStatus status;

		// URL where job status is found.
		string joburl;

		// Last time the jobs status was checked.
		// We rate limit this check.
		chrono::time_point<std::chrono::system_clock> last_poll;
};

// Manages interaction with the CRA.
class CRAHandle {
public:
	// Initializes authentication and CRA application ID.
	CRAHandle() :
		username(getenv("TS_CRA_USERNAME")),
		password(getenv("TS_CRA_PASSWORD")),
		active_jobs(0)
		{};
	bool submit_jobs(string filename);
	bool submit_raxml(CRAJob& job);
	bool parse_status(CRAJob& job);
private:
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

	// Number of CRA jobs currently active.
	// This must stay below some threshold.
	int active_jobs;

	// List of jobs to run.
	vector<CRAJob> jobs;
};
