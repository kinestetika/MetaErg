from time import sleep
from ftplib import FTP, error_perm
from pathlib import Path
from hashlib import md5

NCBI_FTP = None


def fetch(accession, targets: list, checksum_dir: Path) -> dict[Path, bool]:
    # targets is a list of tuples (extension, destination_file)
    outcomes = {destination_file: False for extension, destination_file in targets}
    global NCBI_FTP
    if not NCBI_FTP:
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    local_checksum_file = checksum_dir / 'md5checksums.txt'
    local_checksum_file.unlink(missing_ok=True)
    actual_checksums = {}
    try:
        NCBI_FTP.cwd('/genomes/all/')
        for acc_part in (accession[0:3], accession[4:7], accession[7:10], accession[10:13]):
            NCBI_FTP.cwd(acc_part)
        ftp_dir_list = []
        NCBI_FTP.dir('.', ftp_dir_list.append)
        for directory in ftp_dir_list:
            # print(directory.split()[-1])
            # dr-xr-xr-x   2 ftp      anonymous     4096 Aug 10 05:08 GCF_000295935.1_ASM29593v1
            # dr-xr-xr-x   2 ftp      anonymous     4096 Aug 10 12:42 GCF_000295935.2_IndiAl02
            NCBI_FTP.cwd(directory.split()[-1])
            ftp_file_list = []
            NCBI_FTP.dir('.', ftp_file_list.append)
            targets_found = 0
            checksum_found = 0
            for line in ftp_file_list:
                filename = line.split()[-1]
                for extension, destination_file in targets:
                    if filename.endswith(extension):
                        targets_found += 1
                if 'md5checksums.txt' == filename:
                    checksum_found += 1
            if targets_found and checksum_found:  # the folder contains what we need, proceed...
                for line in ftp_file_list:
                    filename = line.split()[-1]
                    if 'md5checksums.txt' == filename:
                        with open(local_checksum_file, "wb") as local_handle:
                            NCBI_FTP.retrbinary("RETR " + filename, local_handle.write)
                    else:
                        for extension, destination_file in targets:
                            if filename.endswith(extension):
                                with open(destination_file, "wb") as local_handle:
                                    NCBI_FTP.retrbinary("RETR " + filename, local_handle.write)
                                actual_checksums[filename] = {
                                    'md5sum': md5(open(destination_file, 'rb').read()).hexdigest(),
                                    'destination_file': destination_file}
                                outcomes[destination_file] = True
                if sum(outcomes.values()):
                    if local_checksum_file.exists():
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
                                        print(f'Warning: Download failed for {accession}, {filename}: md5sum incorrect.')
                                        Path(actual_checksums[filename]['destination_file']).unlink(missing_ok=True)
                                        outcomes[actual_checksums[filename]['destination_file']] = False
                                except:
                                    pass
                    else:
                        print(f'Warning: could not obtain checksums for {accession}')
            else:
                NCBI_FTP.cwd('..')
        if not sum(outcomes.values()):
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
