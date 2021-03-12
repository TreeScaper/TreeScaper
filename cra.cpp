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
#include "cra.hpp"

using namespace std;

const string cra_url = "https://cipresrest.sdsc.edu/cipresrest/v1/job/";
const string job_completed_string = "COMPLETED";
const string job_submitted_message = "\nYour job has been submitted. You should receive an email at the address "
	"associated with your CIPRES REST API (CRA) account when it is complete. Once finished, re-run TreeScaper, providing the jobname below to perform "
	"tree analysis on the resulting bootstrap trees.\n\n"
	"    ./CLVTreeScaper -trees -jobname ";
const string bootstrap_results_filename = "RAxML_bootstrap.result";

const int poll_time_s = 10;
const int max_poll = 5;

// Populate user buffer with response data
size_t write_callback(char *ptr, size_t size, size_t nmemb, void *userdata) {
	*(string*)userdata += string(ptr);
	return nmemb;
}

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

bool CRAHandle::retrieve_raxml(string jobname) {
	string results_url = cra_url + username + "/" + jobname + "/output";

	CURL *curl = create_cra_curl_handle();

	curl_easy_setopt(curl, CURLOPT_URL, results_url.c_str());

	cout << "Requesting results from job " << jobname << "..." << endl;

	CURLcode res = curl_easy_perform(curl);
	curl_easy_cleanup(curl);

	// Check for errors
	if(res != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));
		return false;
	}

	// Parse response into a pugi xml result
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		cerr << "Error parsing XML" << endl;
		return false;
	}

	for (pugi::xml_node node: doc.child("results").child("jobfiles").children("jobfile")) {
		string filename = node.child("filename").child_value();
		if (filename != bootstrap_results_filename) {
			continue;
		}
		// If filename matches, create a new curl request and download it
		string file_url = node.child("downloadUri").child("url").child_value();
		curl = create_cra_curl_handle();
		curl_easy_setopt(curl, CURLOPT_URL, file_url.c_str());
		res = curl_easy_perform(curl);
		curl_easy_cleanup(curl);

		// Check for errors
		if(res != CURLE_OK) {
			fprintf(stderr, "curl_easy_perform() failed: %s\n",
					curl_easy_strerror(res));
			return false;
		}

		// Write response data to file
		ofstream outfile(bootstrap_results_filename);
		if (outfile) {
			cout << "Writing results to " << bootstrap_results_filename << "..." << endl;
			outfile << userdata;
		}
	}
	return true;
}

bool CRAHandle::submit_raxml(string filename) {

	// Initialize CURL object
	CURL *curl = create_cra_curl_handle();
	if(!curl) {
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
	curl_mime_filedata(field, filename.c_str());

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

	cout << "Submitting job to RAxML on XSEDE..." << endl;

	// Submit job
	CURLcode res = curl_easy_perform(curl);

	curl_easy_cleanup(curl);

	// Check for errors
	if(res != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));
		return false;
	}

	//
	// The code below handles the response XML data
	//

	// Parse response into a pugi xml result
	pugi::xml_document doc;

	pugi::xml_parse_result result = doc.load_buffer(userdata.c_str(), userdata.size());
	if (!result) {
		return false;
	}

	// Get URL to query job status
	string results_url = string(doc.child("jobstatus").child("selfUri").child("url").child_value());

	std::regex r("[^/]*$");
	std::smatch m;
	std::regex_search(results_url, m, r);
	string match = m[0];
	cout << job_submitted_message << match << " <further options>" << endl << endl;

	return true;
}
