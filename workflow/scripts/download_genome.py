import requests
import sys
import os
import gzip
import shutil


sys.stdout = open(snakemake.log[0], 'a')
sys.stderr = sys.stdout


def download_file(url, local_filename):
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
    except requests.exceptions.HTTPError as errh:
        print("HTTP Error: {}".format(errh))
    except requests.exceptions.ConnectionError as errc:
        print("Error Connecting: {}".format(errc))
    except requests.exceptions.Timeout as errt:
        print("Timeout Error: {}".format(errt))
    except requests.exceptions.RequestException as err:
        print("RequestException: {}".format(err))
    except Exception as e:
        print("An error occurred: {}".format(e))



def unzip_file(gz_file):
    unzipped_file = gz_file.replace('.gz', '')
    with gzip.open(gz_file, 'rb') as f_in:
        with open(unzipped_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(gz_file)
    return unzipped_file


assembly = snakemake.params.assembly
provider = snakemake.params.provider
output = snakemake.output[0]

#prepare the download url
base_url = "http://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/{}.fa.gz"
genome_url = base_url.format(assembly, assembly)

#prepare the output file
outFileGZ = output + ".gz"

if not os.path.exists(output):
    print("Downloading {}...".format(genome_url))
    download_file(genome_url, outFileGZ)
    print("Downloaded to {}".format(outFileGZ))
    print("Unzipping {}...".format(outFileGZ))
    unzipped_genome= unzip_file(outFileGZ)
    print("Unzipped to {}".format(unzipped_genome))
else:
    print("{} already exists. Skipping download.".format(output))


