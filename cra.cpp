// Interfaces with CIPRES REST API

#include <string>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <regex>
#include <fstream>
#include <curl/curl.h>
#include <pugixml.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include "cra.hpp"

using namespace std;

// Base URL for CRA API
const string cra_url = "https://cipresrest.sdsc.edu/cipresrest/v1/job/";

// RAxML boostrap results name
const string bootstrap_results_filename = "RAxML_bootstrap.result";

// Directory where results are downloaded
const string cra_output_directory = "cra_output";

// File that tracks status of jobs for each input file
const string job_status_filepath = cra_output_directory + "/job_status";

// Number of jobs that can be active via the CRA at once. This is a limit build
// into the CRA, but we reproduce it here.
const int active_job_limit = 50;

// Minimum time in milliseconds between REST requests
const int rate_limit = 500;

// Minimum time in milliseconds between requests checking status of job.
const int status_rate_limit = 10 * 1000;

// Characters output indicating status of each job in the status file
const map<enum JobStatus, string> status_characters = {
	{UNSUBMITTED, "U"},
	{SUBMITTED, "S"},
	{SUCCESSFUL, "C"},
	{FAILED, "F"},
};

//
// Populate user buffer with response data
//
size_t write_callback(char *ptr, size_t size, size_t nmemb, void *userdata) {
	*(string*)userdata += string(ptr);
	return nmemb;
}

//
// Create a simple CURL handle with the basic configuration needed
// for all CRA requests.
//
CURL *CRAHandle::create_cra_curl_handle() {
	CURL *curl = curl_easy_init();
	if (!curl) {
		return NULL;
	}

	// Set up authentication
	curl_easy_setopt(curl, CURLOPT_USERNAME, username.c_str());
	curl_easy_setopt(curl, CURLOPT_PASSWORD, password.c_str());

	// Create a new list of headers
	struct curl_slist *headerlist = NULL;
	headerlist = curl_slist_append(headerlist, "");

	// Add application key to header list
	string key_header = "cipres-appkey:" + appkey;
	headerlist = curl_slist_append(headerlist, key_header.c_str());

	// Add Expect: 100-continue to header list
	headerlist = curl_slist_append(headerlist, "Expect: 100-continue");

	// Provide header list to curl object
	curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headerlist);

	// Verbose output for debugging. Not desired for production.
	//curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

	// Clear userdata
	userdata = "";

	// Set up callback function and buffer for initial request
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void*) &userdata);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);

	return curl;
}

//
// Performs a request to the CRA via the given URL.
//
// Simple rate-limiting is included, where no more than 1 request per <rate_limit> ms
// may be made.
//
bool CRAHandle::retrieve_url(string url) {

	CURL *curl = create_cra_curl_handle();

	curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

	using namespace std::chrono;

	// Rate limiting
	static time_point<system_clock> last_request_time;

	// Calculate the difference in milliseconds since the last request was made.
	long diff = duration_cast<milliseconds>(system_clock::now() - last_request_time).count();

	// If this difference is less than the limit, sleep for the remainder.
	if (diff < rate_limit) {
		std::this_thread::sleep_for(std::chrono::milliseconds(diff));
	}

	// Perform the request.
	CURLcode res = curl_easy_perform(curl);

	// Reset the last request time.
	last_request_time = system_clock::now();

	curl_easy_cleanup(curl);

	// Check for errors
	if(res != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));
		return false;
	}
	return true;
}

//
// Parses results from the last request to the job url.  These results should already
// be present in the member userdata variable.
//
// If there has been a change in the job's status, change job.status.
// If the job completed successfully, download its output bootstrap trees.
//
// This function expects a specific XML scheme from the CRA.
//
bool CRAHandle::parse_status(CRAJob& job) {

	// Parse response into a pugi xml doc
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		cerr << "Error parsing XML" << endl;
		return false;
	}

	// Check if completed

	bool completed = false;
	bool failed = false;

	// Go through each message and check the "stage". If we find COMPLETED then the
	// job is finished. We then check the <failed> tag for success.
	for (pugi::xml_node node: doc.child("jobstatus").child("messages").children("message")) {
		string stage = node.child("stage").child_value();
		if (stage == "COMPLETED") {
			completed = true;
			string failed_tag =
				string(doc.child("jobstatus").child("failed").child_value());
			if (failed_tag == "true") {
				failed = true;
			}
			break;
		}
	}

	// If completed, first set the job status to either FAILED or SUCCESSFUL.
	// If successful, then download the output file.
	if (!completed) {
		return true;
	}

	// The job has completed, so decrement the active_job counter.
	active_jobs--;

	// If failed, exit before attempting to retrieve results.
	if (failed) {
		change_job_status(job, FAILED);
		return true;
	}

	// If completed and not failed, the job is considered SUCCESSFUL.
	change_job_status(job, SUCCESSFUL);

	// Retrieve list of output files
	string results_url = string(doc.child("jobstatus").child("resultsUri").child("url").child_value());
	retrieve_url(results_url);

	// Parse response into a pugi xml doc
	result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		cerr << "Error parsing XML" << endl;
		return false;
	}

	// The document lists all of the output files, each with a name and URL for retrieving the file.
	// Search through these files until we find the bootstrap results, and then download it with the
	// associated URL.
	for (pugi::xml_node node: doc.child("results").child("jobfiles").children("jobfile")) {

		// Check if name matches expected name for bootstrap results.
		string filename = node.child("filename").child_value();
		if (filename != bootstrap_results_filename) {
			continue;
		}

		// If filename matches, create a new curl request and download it
		string file_url = node.child("downloadUri").child("url").child_value();
		CURL *curl = create_cra_curl_handle();
		curl_easy_setopt(curl, CURLOPT_URL, file_url.c_str());
		CURLcode res = curl_easy_perform(curl);
		curl_easy_cleanup(curl);

		// Check for errors
		if(res != CURLE_OK) {
			fprintf(stderr, "curl_easy_perform() failed: %s\n",
					curl_easy_strerror(res));
			return false;
		}

		// Write response data to file
		// The output files are named based on the input filenames.
		string inputfile_name = job.inputfile.substr(job.inputfile.find_last_of("/\\") + 1);
		string output_filepath = cra_output_directory + "/" + inputfile_name + "_" + bootstrap_results_filename;
		ofstream outfile(output_filepath);
		if (outfile) {
			outfile << userdata;
		}
	}
	return true;
}

//
// Writes a list with the status of each job. Each row is space-delimited with the format:
//
//	 <job status> <input file> <job url>
//
// The statuses are recorded as single characters, as provided by a map<enum JobStatus, string>
//
bool CRAHandle::write_job_status() {
	ofstream of(job_status_filepath);
	if (!of) {
		cout << "Error opening job status path." << endl;
		return false;
	}

	// For each job, write its status, file, and url if it exists
	for (CRAJob& job : jobs) {
		string status_character = status_characters.at(job.status);
		of << status_character << " " << job.inputfile << " " << job.joburl << endl;
	}
	return true;
}

//
// Updates the status variable of a job and writes all job status to the status file.
// This is provided because the status file should be updated every time a jobs status updates.
//
bool CRAHandle::change_job_status(CRAJob& job, enum JobStatus status) {
	job.status = status;
	if (!write_job_status()) {
		return false;
	}
	return true;
}

//
// Handles the main logic for submitting a list of jobs.
//
// After some initialization, it runs a main loop until all jobs have been completed.
// This main loop consists of two inner loops:
//
// The first runs only until the active job limit has not been reached. It goes through all
// jobs submitting unsubmitted jobs.
//
// The second goes through all jobs checking and updating the status for submitted jobs.
//
bool CRAHandle::submit_jobs(string filename) {

	// Create CRAJob objects from input files
	ifstream f(filename);
	string line;
	while (getline(f, line)) {
		jobs.push_back(CRAJob(line));
	}

	// Create output directory if needed
	bool output_directory_exists = false;
	struct stat info;

	// If stat can access the path and it is a directory, then we don't need to recreate
	if (stat(cra_output_directory.c_str(), &info) == 0 && info.st_mode & S_IFDIR) {
		output_directory_exists = true;
	}

	if (!output_directory_exists) {
		// Solution from https://stackoverflow.com/questions/20358455/cross-platform-way-to-make-a-directory-including-subfolders
		mode_t output_directory_mode = 0733;
		int mkdir_res = 0;
		mkdir_res = mkdir(cra_output_directory.c_str(), output_directory_mode);
		if (mkdir_res != 0) {
			cout << "Error creating output directory." << endl;
			return false;
		}
	}

	bool all_jobs_completed = false;

	// Main submission loop
	while (!all_jobs_completed) {

		// Loop through jobs and submit unsubmitted jobs.
		for (CRAJob& job : jobs) {
			if (active_jobs == active_job_limit) {
				break;
			}

			// If we have an unsubmitted job, submit it.
			if (job.status == UNSUBMITTED) {
				submit_raxml(job);
			}
		}

		all_jobs_completed = true;

		using namespace std::chrono;

		// Loop through jobs and check status of submitted jobs.
		for (CRAJob& job : jobs) {
			const time_point<system_clock> now = system_clock::now();
			if (job.status == SUBMITTED &&
				duration_cast<milliseconds>(now - job.last_poll).count() > status_rate_limit) {

				// Get status
				retrieve_url(job.joburl);
				job.last_poll = system_clock::now();
				parse_status(job);
			}

			// Mark that we still have some unfinished jobs.
			if (job.status != SUCCESSFUL && job.status != FAILED) {
				all_jobs_completed = false;
			}
		}
	}

	return true;
}

// Submits a RAxML job
bool CRAHandle::submit_raxml(CRAJob& job) {

	// Initialize CURL object
	CURL *curl = create_cra_curl_handle();
	if (!curl) {
		return false;
	}

	// Set up URL for initial request
	curl_easy_setopt(curl, CURLOPT_URL, (cra_url + username).c_str());

	//
	// The following lines send the initial POST request to submit the CRA job.
	//
	curl_mime *form = NULL;
	curl_mimepart *field = NULL;

	// Create a mime form
	form = curl_mime_init(curl);

	// Fill in the file upload field
	field = curl_mime_addpart(form);
	curl_mime_name(field, "input.infile_");
	curl_mime_filedata(field, job.inputfile.c_str());

	// Specify the CRA tool
	field = curl_mime_addpart(form);
	curl_mime_name(field, "tool");
	curl_mime_data(field, "RAXMLHPC8_REST_XSEDE", CURL_ZERO_TERMINATED);

	// Specify a status email should be sent
	field = curl_mime_addpart(form);
	curl_mime_name(field, "metadata.statusEmail");
	curl_mime_data(field, "true", CURL_ZERO_TERMINATED);

	// Set -x flag for RAxML to perform bootstrap analysis
	field = curl_mime_addpart(form);
	curl_mime_name(field, "vparam.choose_bootstrap_");
	curl_mime_data(field, "b", CURL_ZERO_TERMINATED);

	// Print branch lengths on bootstrap trees
	field = curl_mime_addpart(form);
	curl_mime_name(field, "vparam.printbrlength_");
	curl_mime_data(field, "1", CURL_ZERO_TERMINATED);

	// Provide form to curl object
	curl_easy_setopt(curl, CURLOPT_MIMEPOST, form);

	// Submit job
	CURLcode res = curl_easy_perform(curl);

	curl_easy_cleanup(curl);

	// Check for errors
	if(res != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));
		return false;
	}

	// Job has been submitted successfully, increment number of active_jobs.
	active_jobs++;

	// Parse response into a pugi xml result
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		return false;
	}

	// Get URL to query job status
	string joburl = string(doc.child("jobstatus").child("selfUri").child("url").child_value());
	job.joburl = joburl;

	// Wait until we've updated the joburl to change status, s.t. that status file contains the url.
	change_job_status(job, SUBMITTED);

	return true;
}
