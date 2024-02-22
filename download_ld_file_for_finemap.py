import os
import requests


def download_file(url, path):
    """Downloads a file from the given URL to the specified path."""
    response = requests.get(url)
    if response.status_code == 200:
        with open(path, 'wb') as file:
            file.write(response.content)
    else:
        print(f"Failed to download {url}. Status code: {response.status_code}")


def check_and_download(download_path, chromosome, start, end):
    """Check if the files exist and download them if they don't."""
    file_a_name = f"chr{chromosome}_{start}_{end}.gz"
    file_b_name = f"chr{chromosome}_{start}_{end}.npz"

    file_a_path = os.path.join(download_path, file_a_name)
    file_b_path = os.path.join(download_path, file_b_name)

    base_url = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD"

    # Check and download file A
    if not os.path.exists(file_a_path):
        file_a_url = f"{base_url}/{file_a_name}"
        download_file(file_a_url, file_a_path)
        print(f"Downloaded: {file_a_name}")
    else:
        print(f"File already exists: {file_a_name}")

    # Check and download file B
    if not os.path.exists(file_b_path):
        file_b_url = f"{base_url}/{file_b_name}"
        download_file(file_b_url, file_b_path)
        print(f"Downloaded: {file_b_name}")
    else:
        print(f"File already exists: {file_b_name}")


# Example usage
download_path = "downloadpath"  # Specify the download path
chromosomes = range(1, 23)  # Chromosomes 1 to 22
start_end_pairs = zip(range(1, 250000001, 1000000), range(3000001, 253000001, 1000000))

# Create the download directory if it doesn't exist
if not os.path.exists(download_path):
    os.makedirs(download_path)

# Iterate through chromosomes and start-end pairs
for chromosome in chromosomes:
    for start, end in start_end_pairs:
        check_and_download(download_path, chromosome, start, end)
