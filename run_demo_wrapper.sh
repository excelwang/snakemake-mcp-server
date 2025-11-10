#!/bin/bash
set -x # Enable shell debugging

# This script executes the samtools/faidx wrapper via the /tool-processes REST API endpoint
# using a curl command, polls its status, and verifies successful execution.
# It dynamically fetches the demo case from the /demo-case endpoint.
#
# Prerequisites:
# 1. FastAPI server running on http://127.0.0.1:8082 (e.g., by running 'swa rest --port 8082' in a separate terminal).
# 2. 'jq' command-line JSON processor installed.
# 3. 'curl' installed.

# --- Configuration ---
API_SERVER_URL="http://127.0.0.1:8082"
MAX_ATTEMPTS=60 # 60 seconds timeout for job polling

# --- Step 1: Fetch demo case from /demo-case endpoint ---
echo "--- Step 1: Fetching demo case from /demo-case endpoint ---"
DEMO_CASE_RESPONSE=$(curl -s "$API_SERVER_URL/demo-case")

if [ $? -ne 0 ]; then
    echo "Error: Failed to connect to API server at $API_SERVER_URL. Is it running?"
    exit 1
fi

DEMO_CASE_METHOD=$(echo "$DEMO_CASE_RESPONSE" | jq -r '.method')
DEMO_CASE_ENDPOINT=$(echo "$DEMO_CASE_RESPONSE" | jq -r '.endpoint')
PAYLOAD_JSON=$(echo "$DEMO_CASE_RESPONSE" | jq -r '.payload | tojson')
CURL_EXAMPLE=$(echo "$DEMO_CASE_RESPONSE" | jq -r '.curl_example')

if [ "$DEMO_CASE_METHOD" != "POST" ] || [ "$DEMO_CASE_ENDPOINT" != "/tool-processes" ]; then
    echo "Error: Unexpected demo case method or endpoint."
    echo "$DEMO_CASE_RESPONSE" | jq .
    exit 1
fi

echo "Fetched Demo Case:"
echo "$DEMO_CASE_RESPONSE" | jq .
echo ""

# Extract details from the payload
WORKDIR=$(echo "$PAYLOAD_JSON" | jq -r '.workdir')
INPUT_FILE_RELATIVE=$(echo "$PAYLOAD_JSON" | jq -r '.inputs[0]')
OUTPUT_FILE_RELATIVE=$(echo "$PAYLOAD_JSON" | jq -r '.outputs[0]')

echo "Extracted Workdir: $WORKDIR"
echo "Extracted Input File: $INPUT_FILE_RELATIVE"
echo "Extracted Output File: $OUTPUT_FILE_RELATIVE"
echo ""

# --- Step 2: Set up environment and create dummy input file ---
echo "--- Step 2: Setting up environment and creating dummy input file ---"

# Ensure the workdir exists (it should be created by the /demo-case endpoint, but good to be safe)
mkdir -p "$WORKDIR"

# Create a dummy input file in the demo directory with content matching test_direct_fastapi_demo_case_endpoint
echo ">chr1" > "$WORKDIR/$INPUT_FILE_RELATIVE"
echo "AGCTAGCTAGCTAGCT" >> "$WORKDIR/$INPUT_FILE_RELATIVE"
echo ">chr2" >> "$WORKDIR/$INPUT_FILE_RELATIVE"
echo "TCGATCGATCGA" >> "$WORKDIR/$INPUT_FILE_RELATIVE"

echo "Created input file: $WORKDIR/$INPUT_FILE_RELATIVE"
echo "Expected output file: $WORKDIR/$OUTPUT_FILE_RELATIVE"
echo ""

# --- Step 3: Call the /tool-processes endpoint to start the wrapper execution ---
echo "--- Step 3: Calling /tool-processes endpoint ---"
read -r CURL_OUTPUT < <(curl -s -X POST "$API_SERVER_URL$DEMO_CASE_ENDPOINT" \
     -H "Content-Type: application/json" \
     -d "$PAYLOAD_JSON")

echo "Response from /tool-processes:"
echo "$CURL_OUTPUT" | jq .
echo ""

# --- Step 4: Extract job_id and status_url from the response ---
echo "--- Step 4: Extracting job_id and status_url ---"
JOB_ID=$(echo "$CURL_OUTPUT" | jq -r '.job_id')
STATUS_URL_RELATIVE=$(echo "$CURL_OUTPUT" | jq -r '.status_url')
STATUS_URL="$API_SERVER_URL$STATUS_URL_RELATIVE"

if [ -z "$JOB_ID" ] || [ "$JOB_ID" == "null" ]; then
    echo "Error: Failed to get JOB_ID from response."
    echo "Full response: $CURL_OUTPUT"
    exit 1
fi

echo "Job ID: $JOB_ID"
echo "Status URL: $STATUS_URL"
echo ""

# --- Step 5: Poll the status_url to monitor job progress ---
echo "--- Step 5: Polling job status ---"
JOB_STATUS=""
ATTEMPT=0

while [ "$JOB_STATUS" != "completed" ] && [ "$JOB_STATUS" != "failed" ] && [ "$ATTEMPT" -lt "$MAX_ATTEMPTS" ]; do
    echo "Polling job status... Attempt $((ATTEMPT+1)) of $MAX_ATTEMPTS"
    sleep 1
    JOB_RESPONSE=$(curl -s "$STATUS_URL")
    JOB_STATUS=$(echo "$JOB_RESPONSE" | jq -r '.status')
    echo "Current job status: "$JOB_STATUS""
    ATTEMPT=$((ATTEMPT+1))
done

echo ""
echo "Final job response:"
echo "$JOB_RESPONSE" | jq .
echo ""

# --- Step 6: Verify the final job status and output file ---
echo "--- Step 6: Verifying final job status and output ---"
if [ "$JOB_STATUS" == "completed" ]; then
    echo "SUCCESS: Wrapper execution completed successfully!"
    # Add a small delay to ensure file system reflects the change
    sleep 2
    # Check if the output file exists
    if [ -f "$WORKDIR/$OUTPUT_FILE_RELATIVE" ]; then
        echo "Output file $WORKDIR/$OUTPUT_FILE_RELATIVE created successfully."
        echo "Content of output file:"
        cat "$WORKDIR/$OUTPUT_FILE_RELATIVE"
        
        # Verify content matches the expected output from test_direct_fastapi_demo_case_endpoint
        EXPECTED_CONTENT="chr1\t16\t6\t16\t17\nchr2\t12\t29\t12\t13\n"
        ACTUAL_CONTENT=$(cat "$WORKDIR/$OUTPUT_FILE_RELATIVE")

        if [ "$ACTUAL_CONTENT" == "$EXPECTED_CONTENT" ]; then
            echo "Output file content verified successfully."
            exit 0
        else
            echo "ERROR: Output file content mismatch."
            echo "Expected:"
            echo -e "$EXPECTED_CONTENT"
            echo "Actual:"
            echo -e "$ACTUAL_CONTENT"
            exit 1
        fi
    else
        echo "ERROR: Output file $WORKDIR/$OUTPUT_FILE_RELATIVE was not created."
        exit 1
    fi
else
    echo "ERROR: Wrapper execution failed or timed out."
    echo "Job details:"
    echo "$JOB_RESPONSE" | jq '.'
    exit 1
fi