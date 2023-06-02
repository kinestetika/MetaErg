from time import sleep
from ftplib import FTP, error_perm
from pathlib import Path
from hashlib import md5

NCBI_FTP = None

def fetch(accession, targets: list, checksum_dir: Path) -> list[bool]:
    # targets is a list of tuples (extension, destination_file)
    global NCBI_FTP
    if not NCBI_FTP:
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()

    outcomes = []
    local_checksum_file = checksum_dir / 'md5checksums.txt'
    local_checksum_file.unlink(missing_ok=True)
    actual_checksums = {}
    try:
        NCBI_FTP.cwd('/genomes/all/')
        for acc_part in (accession[0:3], accession[4:7], accession[7:10], accession[10:13]):
            NCBI_FTP.cwd(acc_part)
        ftp_dir_list = []
        NCBI_FTP.dir('.', ftp_dir_list.append)
        NCBI_FTP.cwd(ftp_dir_list[0].split()[-1])
        ftp_dir_list = []
        NCBI_FTP.dir('.', ftp_dir_list.append)
        # print('actual_checksums:')
        for extension, destination_file in targets:
            for l in ftp_dir_list:
                filename = l.split()[-1]
                if filename.endswith(extension):
                    with open(destination_file, "wb") as local_handle:
                        NCBI_FTP.retrbinary("RETR " + filename, local_handle.write)
                    actual_checksums[filename] = {'md5sum': md5(open(destination_file, 'rb').read()).hexdigest(),
                                                  'destination_file': destination_file}
                    # print(actual_checksums[filename], filename)
                    outcomes.append(filename)
                elif 'md5checksums.txt' == filename:
                    with open(local_checksum_file, "wb") as local_handle:
                        NCBI_FTP.retrbinary("RETR " + filename, local_handle.write)
        if outcomes:
            if local_checksum_file.exists():
                # print('expected_checksums:')
                with open(local_checksum_file) as checksum_content:
                    for line in checksum_content:
                        expected_checksum, filename = line.split()
                        filename = filename[2:]
                        try:
                            actual_checksum = actual_checksums[filename]['md5sum']
                            if actual_checksum == expected_checksum:
                                # print(expected_checksum, filename, 'OK')
                                pass
                            else:
                                print(f'Warning: Download failed for {accession}: md5sum incorrect.')
                                Path(actual_checksums[filename]['destination_file']).unlink(missing_ok=True)
                                outcomes = []
                        except:
                            pass
            else:
                print(f'Warning: could not obtain checksums for {accession}')
        else:
            print(f'Warning: no annotations at NCBI for {accession}')
        return outcomes
    except EOFError:
        print('FTP Error at NCBI - resetting connection')
        sleep(10)
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    except BrokenPipeError:
        print('FTP Error at NCBI - resetting connection')
        sleep(10)
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    except ConnectionRefusedError:
        print('FTP Error at NCBI - resetting connection')
        sleep(10)
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    except error_perm:
        print(f'Unable to find {accession}')
        return outcomes

def main():
    fetch('GCA_020633295', [('_protein.faa.gz', 'protein_file.faa.gz'), ('_rna_from_genomic.fna.gz', 'rna_file.faa.gz')], Path())

if __name__ == "__main__":
    main()
