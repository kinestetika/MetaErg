from ftplib import FTP, error_perm

NCBI_FTP = None


def fetch(accession, targets: list) -> list[bool]:
    # targets is a list of tuples (extension, destination_file)
    global NCBI_FTP
    if not NCBI_FTP:
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()

    outcomes = []
    try:
        NCBI_FTP.cwd('/genomes/all/')
        for acc_part in (accession[0:3], accession[4:7], accession[7:10], accession[10:13]):
            NCBI_FTP.cwd(acc_part)
        ftp_dir_list = []
        NCBI_FTP.dir('.', ftp_dir_list.append)
        NCBI_FTP.cwd(ftp_dir_list[0].split()[-1])
        ftp_dir_list = []
        NCBI_FTP.dir('.', ftp_dir_list.append)
        for extension, destination_file in targets:
            success = False
            for l in ftp_dir_list:
                filename = l.split()[-1]
                if filename.endswith(extension):
                    success = True
                    with open(destination_file, "wb") as local_handle:
                        NCBI_FTP.retrbinary("RETR " + filename, local_handle.write)
            outcomes.append(success)
        return outcomes
    except EOFError:
        print('FTP Error at NCBI - resetting connection')
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    except BrokenPipeError:
        print('FTP Error at NCBI - resetting connection')
        NCBI_FTP = FTP('ftp.ncbi.nlm.nih.gov')
        NCBI_FTP.login()
    except error_perm:
        print(f'Unable to find {accession}')
        return outcomes