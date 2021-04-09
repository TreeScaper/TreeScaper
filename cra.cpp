// Interfaces with CIPRES REST API

#include <string>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <regex>
#include <fstream>
#include <iomanip>
#include <curl/curl.h>
#include <pugixml.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include "cra.hpp"

using namespace std;

// Level for logging to stdout
enum CRALogLevel cra_log_level;

// Base URL for CRA API
const string cra_url = "https://cipresrest.sdsc.edu/cipresrest/v1/job/";

// Parameter for name of the output file to download.
const string output_file_param = "treescaper.output";

// Directory where results are downloaded
const string cra_output_directory = "cra_output";

// File that tracks status of jobs for each input file
const string job_status_filepath = cra_output_directory + "/job_status";

// Delimiter between fields in status file.
const char status_file_delim = '\t';

// Number of jobs that can be active via the CRA at once. This is a limit build
// into the CRA, but we reproduce it here.
const int active_job_limit = 50;

// Minimum time in milliseconds between REST requests
const int rate_limit = 2000;

// Application ID for CRA
const string cra_application_id = "treescaper_inference_dev-D4DFA6E180C643779DC203C1D5114ED4";

// Increase the poll interval from the minimum allowable as a matter of courtesy.
const int poll_interval_multiplier = 2;

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
// Populate buffer with header data
//
size_t header_callback(char *ptr, size_t size, size_t nmemb, void *headerdata) {
	*(string*)headerdata += string(ptr);
	return nmemb;
}

void log_message(string message) {
	if (cra_log_level == DEBUG) {
		time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		cout << std::put_time(gmtime(&time), "CRA %Y-%m-%d %H:%M:%S") << ": " << message << endl;
	}
}

//
// Submits request and checks for errors and bad response codes.
//
bool submit_curl_request(CURL *curl) {

	using namespace std::chrono;

	// Rate limiting
	static time_point<system_clock> last_request_time;

	// Calculate the difference in milliseconds since the last request was made.
	long diff = duration_cast<milliseconds>(system_clock::now() - last_request_time).count();

	// If this difference is less than the limit, sleep for the remainder.
	if (diff < rate_limit) {
		std::this_thread::sleep_for(std::chrono::milliseconds(rate_limit - diff));
	}

	log_message("Sending request..");

	// Perform the request.
	CURLcode res = curl_easy_perform(curl);

	curl_easy_cleanup(curl);

	// Check for errors
	if(res != CURLE_OK) {
		cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << endl;
		return false;
	}

	long response_code = 0;
	curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);

	// Check for HTTP OK response code.
	if (response_code != 200) {
		log_message("FAIL");
		cerr << "Error: non-200 HTTP response. Response code was " << response_code << "." << endl;
		return false;
	}
	log_message("OK");

	// Reset the last request time.
	last_request_time = system_clock::now();

	return true;
}

//
// Logs data and headers from most recent HTTP response.
//
void CRAHandle::log_last_response() {
	if (cra_log_level = DEBUG) {
		cout << "Response data:" << endl;
		cout << userdata << endl;
		cout << "Headers:" << endl;
		cout << headerdata << endl;
	}
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
	string key_header = "cipres-appkey:" + cra_application_id;
	headerlist = curl_slist_append(headerlist, key_header.c_str());

	// Add Expect: 100-continue to header list
	headerlist = curl_slist_append(headerlist, "Expect: 100-continue");

	// Provide header list to curl object
	curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headerlist);

	// Verbose output for debugging. Not desired for production.
	//curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

	// Clear response and header data.
	userdata = "";
	headerdata = "";

	// Set up callback function and buffer for initial request
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void*) &userdata);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);

	// Set up callback function and buffer for header data
	curl_easy_setopt(curl, CURLOPT_HEADERDATA, (void*) &headerdata);
	curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, header_callback);

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

	log_message("Setting request URL: " + url);

	curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

	// Submit request.
	if (!submit_curl_request(curl)) {
		log_last_response();
		return false;
	}

	return true;
}

bool CRAHandle::parse_single_status(pugi::xml_node status_node) {

	// Check if job finished.
	bool finished = (string(status_node.child("terminalStage").child_value()) == string("true"));
	if (!finished) {
		return true;
	}

	// Find matching job from list.
	string job_handle = string(status_node.child("jobHandle").child_value());
	CRAJob *job;
	for (CRAJob& job_p : jobs) {
		if (job_p.status == SUBMITTED) {
			if (job_p.handle == job_handle) {
				job = &job_p;
				break;
			}
		}
	}

	// Decrement the number of jobs if the job has finished.
	active_jobs--;

	// Check if there was a failure by CRA to submit the job or retrieve results.
	bool failed = (string(status_node.child("failed").child_value()) == string("true"));
	if (failed) {
		change_job_status(*job, FAILED);
		return false;
	}

	// Check if we have a message that the job completed successfully.
	bool completed_message = false;
	for (pugi::xml_node node: status_node.child("messages").children("message")) {
		string stage = node.child("stage").child_value();
		if (stage == "COMPLETED") {
			completed_message = true;
			break;
		}
	}

	// If we are at the terminal stage but don't have that message, consider it failed.
	if (!completed_message) {
		change_job_status(*job, FAILED);
		return true;
	}

	// Retrieve list of output files
	string results_url = string(status_node.child("resultsUri").child("url").child_value());
	retrieve_url(results_url);

	// Parse response into a pugi xml doc
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		cerr << "Error parsing XML: " << userdata << endl;
		return false;
	}

	// The document lists all of the output files, each with a name and URL for retrieving the file.
	// Search through these files until we find the bootstrap results, and then download it with the
	// associated URL.
	for (pugi::xml_node node: doc.child("results").child("jobfiles").children("jobfile")) {

		// Check if name matches expected name for bootstrap results.
		string filename = node.child("filename").child_value();
		string output_filename = cra_params.at(output_file_param);
		if (filename != output_filename) {
			continue;
		}

		// If filename matches, create a new curl request and download it
		string file_url = node.child("downloadUri").child("url").child_value();
		CURL *curl = create_cra_curl_handle();

		log_message("Setting request URL: " + file_url);

		curl_easy_setopt(curl, CURLOPT_URL, file_url.c_str());

		// Submit request.
		if (!submit_curl_request(curl)) {
			log_last_response();
			return false;
		}

		// Write response data to file
		// The output files are named based on the input filenames.
		string inputfile_name = job->inputfile.substr(job->inputfile.find_last_of("/\\") + 1);
		string output_filepath = cra_output_directory + "/" + inputfile_name + "_" + output_filename;
		ofstream outfile(output_filepath);
		if (outfile) {
			outfile << userdata;
		}
		outfile.close();
	}

	// If completed and not failed, the job is considered SUCCESSFUL.
	change_job_status(*job, SUCCESSFUL);
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
bool CRAHandle::parse_status_list() {

	// Parse response into a pugi xml doc
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		cerr << "Error parsing XML: " << userdata << endl;
		return false;
	}

	// Parse the status of each job
	for (pugi::xml_node node: doc.child("joblist").child("jobs").children("jobstatus")) {
		if (!parse_single_status(node)) {
			return false;
		}
	}

	return true;
}

//
// Writes a list with the status of each job. Each row is tab-delimited with the format:
//
//	 <input file> <job status> <job url>
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
		of << job.inputfile << status_file_delim << status_character << status_file_delim << job.handle << endl;
	}
	of.close();
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
// jobs submitting unsubmitted jobs. The second goes through all jobs checking and updating the status for submitted jobs.
//
// @param filelist  File containing a list of input files.
// @param paramfile File containing common parameters for all jobs.
//
bool CRAHandle::submit_jobs(string filelist, string paramfile) {

	// Parse parameters.
	ifstream param_f(paramfile);
	string param_line;
	while (getline(param_f, param_line)) {

		// Ignore commented and blank lines.
		std::regex comment_regex("^#");
		if (regex_search(param_line, comment_regex) || param_line.length() == 0) {
			continue;
		}

		// Find position of "=".
		size_t delim_pos = param_line.find("=");

		// Parse name and value of configuration parameter
		string param_name = param_line.substr(0, delim_pos);
		string param_value = param_line.substr(delim_pos + 1);

		// Add parameter to param map
		cra_params.insert(std::pair<string,string>(param_name, param_value));
	}
	param_f.close();

	// Create CRAJob objects from input files
	ifstream list_f(filelist);
	string list_line;

	while (getline(list_f, list_line)) {
		istringstream line_stream(list_line);

		// Read filename.
		string inputfile("");
		if (!getline(line_stream, inputfile, status_file_delim) && line_stream.bad()) {
			return false;
		}

		// Read status character.
		string status_character("");
		if (!getline(line_stream, status_character, status_file_delim) && line_stream.bad()) {
			return false;
		}

		// Read URL for job if submitted.
		string job_handle("");
		if (!getline(line_stream, job_handle, status_file_delim) && line_stream.bad()) {
			return false;
		}

		// Convert status character to status enum.
		enum JobStatus status;
		if (status_character == "") {
			status = UNSUBMITTED;
		} else {
			for (pair <enum JobStatus,string> status_pair : status_characters) {
				if (status_pair.second == status_character) {
					status = status_pair.first;
					if (status == SUBMITTED) {
						active_jobs++;
					}
					break;
				}
			}
		}

		CRAJob job(inputfile, status, job_handle);
		jobs.push_back(job);
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

	write_job_status();

	bool all_jobs_completed = false;

	// Main submission loop
	while (!all_jobs_completed) {

		// Loop through jobs and submit unsubmitted jobs.
		for (CRAJob& job : jobs) {
			if (active_jobs >= active_job_limit) {
				break;
			}

			// If we have an unsubmitted job, submit it.
			if (job.status == UNSUBMITTED) {
				if (!submit_job(job)) {
					return false;
				}
			}
		}

		all_jobs_completed = true;

		using namespace std::chrono;

		// Build url to request status of all active jobs
		string url_query = "?";
		int jobs_queried = 0;

		for (CRAJob& job : jobs) {
			if (job.status == SUBMITTED) {
				// The URL cannot be too long, so limit the number of queried jobs to the limit
				// for active jobs. We may have more jobs than this in total to query, because of
				// completed but not active jobs.
				if (jobs_queried++ == active_job_limit){
					break;
				}
				url_query.append("jh=" + job.handle + "&");
			}
		}

		// If enough time has passed, check the status of all submitted jobs.
		const time_point<system_clock> now = system_clock::now();
		if (duration_cast<seconds>(now - last_poll).count() > min_poll_interval_seconds) {

			// Get status of all active jobs
			retrieve_url(cra_url + "/" + url_query);
			last_poll = system_clock::now();

			// Parse status.
			parse_status_list();
		}

		// Check if there are any unfinished jobs.
		for (CRAJob& job : jobs) {
			if (job.status != SUCCESSFUL && job.status != FAILED) {
				all_jobs_completed = false;
				break;
			}
		}
	}

	return true;
}

bool CRAHandle::submit_job(CRAJob& job) {
	// Initialize CURL object
	CURL *curl = create_cra_curl_handle();
	if (!curl) {
		return false;
	}

	log_message("Setting request URL: " + cra_url + username);

	// Set up URL for initial request
	curl_easy_setopt(curl, CURLOPT_URL, (cra_url + username).c_str());

	//
	// The following lines send the initial POST request to submit the CRA job.
	//
	curl_mime *form = NULL;
	curl_mimepart *field = NULL;

	// Create a mime form.
	form = curl_mime_init(curl);

	// Add the input.infile_ parameter separately because it is the only
	// parameter that differs for each job in a batch.
	field = curl_mime_addpart(form);
	curl_mime_name(field, "input.infile_");
	curl_mime_filedata(field, job.inputfile.c_str());

	// Add each parameter to the request.
	for (pair <string,string> param : cra_params) {
		string param_name = param.first;
		string param_value = param.second;

		// Create new field and populate name.
		// The function used to populate the data depends on the type of field.

		std::regex input_regex("^input");
		std::regex treescaper_regex("^treescaper");
		std::regex other_param_regex("^(vparam|tool|metadata)");

		// Parameters prefixed with input expect a file to be uploaded.
		if (regex_search(param_name, input_regex)) {

			// Add filedata field
			field = curl_mime_addpart(form);
			curl_mime_name(field, param_name.c_str());
			curl_mime_filedata(field, param_value.c_str());

		// Parameters prefixed with tool, vparam, or metadata are simply passed as the parameter value.
		} else if (regex_search(param_name, other_param_regex)) {

			// Add data field
			field = curl_mime_addpart(form);
			curl_mime_name(field, param_name.c_str());
			curl_mime_data(field, param_value.c_str(), CURL_ZERO_TERMINATED);

		// TreeScaper-specific configuration is not passed to CRA
		} else if (regex_search(param_name, treescaper_regex)) {
			continue;

		// One of the above prefixes is expected for CRA parameters.
		} else {
			cerr << "Warning: Unknown parameter type for " << param.first << ". Ignoring." << endl;
		}
	}

	// Provide form to curl object
	curl_easy_setopt(curl, CURLOPT_MIMEPOST, form);

	log_message("Submitting job for " + job.inputfile + ".");

	// Submit request.
	if (!submit_curl_request(curl)) {
		log_last_response();
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
	job.handle = string(doc.child("jobstatus").child("jobHandle").child_value());

	log_message("Job handle is: " + job.handle);

	// Get minimum poll interval.
	min_poll_interval_seconds =
		poll_interval_multiplier * doc.child("jobstatus").child("minPollIntervalSeconds").text().as_int();

	// Wait until we've updated the handle to change status, s.t. that status file contains the url.
	change_job_status(job, SUBMITTED);

	return true;
}
