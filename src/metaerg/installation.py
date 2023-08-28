from pathlib import Path
import os
from shutil import which


def install_all_helper_programs(bin_dir: Path, path_to_signalp: Path, path_to_tmhmm: Path, path_to_antismash_db: Path):
    # check for required programs
    success = True
    for cmd in 'git make gcc tar wget perl sed'.split():
        if cmd_loc := which(cmd):
            print(f'{cmd}: will use {cmd_loc}')
        else:
            print(f'{cmd}: not found. Please install first using your system\'s package manager.')
    if not success:
        print('aborted installation of helper programs.')
        return

    # Create the profile file.
    # To "activate" your metaerg installation, you will need to run, for example:
    # >source /home/my_name/bin/metaerg/bin/profile
    # (if that is the path to your installation)
    profile = f'''
export BIOINF_PREFIX={bin_dir}
PATH=$BIOINF_PREFIX/infernal/binaries:$PATH
PATH=$BIOINF_PREFIX/hmmer2/src:$PATH
PATH=$BIOINF_PREFIX/tmhmm/bin:$PATH
PATH=$BIOINF_PREFIX/repeatscout:$PATH
PATH=$BIOINF_PREFIX/repeatmasker:$PATH
PATH=$BIOINF_PREFIX/minced:$PATH
PATH=$BIOINF_PREFIX/hmmer3/bin:$PATH
PATH=$BIOINF_PREFIX/ncbi-blast/bin:$PATH
PATH=$BIOINF_PREFIX/emboss/bin:$PATH
PATH=$BIOINF_PREFIX/vienna_rna/bin:$PATH
PATH=$BIOINF_PREFIX/cd-hit:$PATH
PATH=$BIOINF_PREFIX:$PATH
export PATH
'''
    profile_file = bin_dir / 'profile'
    with open(profile_file, "w") as profile_handle:
        profile_handle.write(profile)
    os.system(f'cat {profile_file}')
    for line in profile.split():
        if line.startswith('PATH='):
            os.environ["PATH"] += ':' + line.replace('PATH=$PATH:$BIOINF_PREFIX', str(bin_dir))
    print(os.environ["PATH"])
    os.chdir(bin_dir)
    # # (minced) minced 0.4.2 https://github.com/ctSkennerton/minced
    # os.system('wget -q https://github.com/ctSkennerton/minced/archive/refs/tags/0.4.2.tar.gz')
    # os.system('tar -xf 0.4.2.tar.gz')
    # os.chdir('minced-0.4.2')
    # os.system('make')
    # os.chdir(bin_dir)
    # os.system('mv minced-0.4.2 minced')
    # os.system('rm 0.4.2.tar.gz')

    # (crisprdetect) crisprdetect 2.4 https://github.com/davidchyou/CRISPRDetect_2.4/tree/master
    # first we need emboss
    os.system("wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz'")
    os.system('mv emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz .')
    os.system('rm -r emboss.open-bio.org')
    os.system('tar -zxf emboss-latest.tar.gz')
    os.system('rm emboss-latest.tar.gz')
    os.chdir('EMBOSS-6.6.0/')
    os.system(f'./configure --prefix={bin_dir / "emboss"}')
    os.system('make')
    os.system(f'make install')
    os.chdir(bin_dir)
    os.system('rm -r EMBOSS-6.6.0/')
    # we also need the vienna RNA suite
    os.system('wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.3.tar.gz')
    os.system('tar -zxf ViennaRNA-2.6.3.tar.gz')
    os.system('rm ViennaRNA-2.6.3.tar.gz')
    os.chdir('ViennaRNA-2.6.3/')
    os.system('./configure --prefix=/home/kinestetika/bin/vienna_rna --without-perl --without-python --without-rnaxplorer')
    os.system('make')
    os.system(f'make install')
    os.chdir(bin_dir)
    os.system('rm -rf ViennaRNA-2.6.3')
    # we also need (cdhit) cd-hit 4.8.1 https://github.com/weizhongli/cdhit
    os.system('wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz')
    os.system('tar -xf cd-hit-v4.8.1-2019-0228.tar.gz')
    os.chdir('cd-hit-v4.8.1-2019-0228')
    os.system('make')
    os.chdir(bin_dir)
    os.system('mv cd-hit-v4.8.1-2019-0228 cd-hit')
    os.system('rm cd-hit-v4.8.1-2019-0228.tar.gz')
    # finally: CRISPRDetect
    os.system('git clone https://github.com/davidchyou/CRISPRDetect_2.4.git')
    os.chdir('CRISPRDetect_2.4/')
    os.system('unzip CRISPRDetect_2.4.zip')
    os.system('rm CRISPRDetect_2.4.zip')
    os.chdir(bin_dir)
    os.system('mv CRISPRDetect_2.4/CRISPRDetect_2.4 CRISPRDetect')
    os.system('rm -rf CRISPRDetect_2.4/')
    os.system('mv CRISPRDetect/clustalw .')
    os.system('chmod a+x clustalw')
    os.system('chmod a+x CRISPRDetect/CRISPRDetect.pl')
    # still need to make sure we install package perl-parallel-forkmanager with pacman

    # (padloc)
    # (aragorn) aragorn 1.2.41 https://www.ansikte.se/ARAGORN/Downloads/
    os.system('wget -q https://www.ansikte.se/ARAGORN/Downloads/aragorn1.2.41.c')
    os.system('gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.41.c')
    os.system('rm aragorn1.2.41.c')
    # (infernal) cmsearch 1.1.4 http://eddylab.org/infernal/
    os.system('wget -q http://eddylab.org/infernal/infernal-1.1.4-linux-intel-gcc.tar.gz')
    os.system('tar -xf infernal-1.1.4-linux-intel-gcc.tar.gz')
    os.system('mv infernal-1.1.4-linux-intel-gcc infernal')
    os.system('rm infernal-1.1.4-linux-intel-gcc.tar.gz')
    # (genometools) gt 1.6.2 http://genometools.org/tools/gt_ltrharvest.html
    os.system('wget -q http://genometools.org/pub/binary_distributions/gt-1.6.2-Linux_x86_64-64bit-barebone.tar.gz')
    os.system('tar -xf gt-1.6.2-Linux_x86_64-64bit-barebone.tar.gz')
    os.system('mv gt-1.6.2-Linux_x86_64-64bit-barebone genometools')
    os.system('ln -sf genometools/bin/gt .')
    os.system('rm gt-1.6.2-Linux_x86_64-64bit-barebone.tar.gz')
    # (tandem-repeat-finder) trf 4.09 https://github.com/Benson-Genomics-Lab/TRF
    os.system('wget -q https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64')
    os.system('chmod a+x trf409.linux64')
    os.system('ln -sf trf409.linux64 trf')
    # (RepeatScout) RepeatScout 1.0.5 https://github.com/mmcco/RepeatScout
    os.system('git clone https://github.com/mmcco/RepeatScout.git')
    os.system('mv RepeatScout repeatscout')
    os.chdir("repeatscout")
    os.system('make')
    os.chdir(bin_dir)
    #(nseg) nseg 1.0.1 https://github.com/jebrosen/nseg
    #(nseg is needed by repeatscout/repeatmasker)
    os.chdir(bin_dir)
    os.system('wget https://github.com/jebrosen/nseg/archive/refs/tags/v1.0.1.tar.gz')
    os.system('tar -xf v1.0.1.tar.gz')
    os.chdir("nseg-1.0.1")
    os.system('make')
    os.chdir(bin_dir)
    os.system('ln -sf nseg-1.0.1/nseg nseg')
    os.system('ln -sf nseg-1.0.1/nmerge nmerge')
    os.system('rm v1.0.1.tar.gz')
    # (rmblast) rmblastn 2.14.0 https://www.repeatmasker.org/rmblast/
    # (rmblast is needed by repeatmasker)
    os.system('wget -q https://www.repeatmasker.org/rmblast/rmblast-2.14.0+-x64-linux.tar.gz')
    os.system('tar -xf rmblast-2.14.0+-x64-linux.tar.gz')
    os.system('mv rmblast-2.14.0 rmblast')
    os.system('rm rmblast-2.14.0+-x64-linux.tar.gz')
    # (RepeatMasker) RepeatMasker 4.1.5 http://www.repeatmasker.org/RepeatMasker/
    os.system('wget -q http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz')
    os.system('tar -xf RepeatMasker-4.1.5.tar.gz')
    os.system('mv RepeatMasker repeatmasker')
    os.chdir("repeatmasker")
    os.system(f'perl ./configure -default_search_engine rmblast -libdir {bin_dir / "repeatmasker" / "Libraries"} '
              f'-rmblast_dir {bin_dir / "rmblast" / "bin"} -trf_prgm {bin_dir / "trf"}')
    os.chdir(bin_dir)
    os.system('rm RepeatMasker-4.1.5.tar.gz')
    # (prodigal) prodigal 2.6.3 https://github.com/hyattpd/Prodigal
    os.system('wget -q https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux')
    os.system('chmod a+x prodigal.linux')
    os.system('ln -sf prodigal.linux prodigal')
    # (diamond) diamond 2.0.15 https://github.com/bbuchfink/diamond
    # (later versions do not work with antismash)
    os.system('wget -q https://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz')
    os.system('tar -xf diamond-linux64.tar.gz')
    os.system('rm diamond-linux64.tar.gz')
    # (ncbi-blast) blastn 2.14.0 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    os.system('wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz')
    os.system('tar -xf ncbi-blast-2.14.0+-x64-linux.tar.gz')
    os.system('mv ncbi-blast-2.14.0+ ncbi-blast')
    os.system('rm ncbi-blast-2.14.0+-x64-linux.tar.gz')
    # (hmmer-3) hmmsearch 3.3.2 http://hmmer.org
    os.system('wget -q http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz')
    os.system('tar -xf hmmer-3.3.2.tar.gz')
    os.chdir('hmmer-3.3.2')
    os.system(f'./configure --prefix={bin_dir / "hmmer-3.3.2"}')
    os.system('make')
    os.system('make install')
    os.chdir('easel')
    os.system('make install')
    os.chdir(bin_dir)
    os.system('ln -s hmmer-3.3.2 hmmer3')
    # (hmmer-2) hmmsearch2 2.3.2 http://hmmer.org/
    # (required by antismash)
    os.system('wget -q http://eddylab.org/software/hmmer/hmmer-2.3.2.tar.gz')
    os.system('tar -xf hmmer-2.3.2.tar.gz')
    os.chdir('hmmer-2.3.2')
    os.system('./configure')
    os.system('make')
    os.chdir('src')
    os.system('mv hmmalign hmmalign2')
    os.system('mv hmmbuild hmmbuild2')
    os.system('mv hmmcalibrate hmmcalibrate2')
    os.system('mv hmmconvert hmmconvert2')
    os.system('mv hmmemit hmmemit2')
    os.system('mv hmmfetch hmmfetch2')
    os.system('mv hmmindex hmmindex2')
    os.system('mv hmmpfam hmmpfam2')
    os.system('mv hmmsearch hmmsearch2')
    os.chdir(bin_dir)
    os.system('ln -s hmmer-2.3.2 hmmer2')
    os.system('rm hmmer-2.3.2.tar.gz')
    #(muscle) muscle 5.1 https://drive5.com/muscle
    #(no longher required by antismash, version 7 and up)
    os.chdir(bin_dir)
    os.system('wget -q https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz')
    os.system('tar -xf muscle3.8.31_i86linux64.tar.gz')
    os.system('ln -sf muscle3.8.31_i86linux64 muscle')
    os.system('rm muscle3.8.31_i86linux64.tar.gz')
    #(fasttree) fasttree 2.1.11 https://microbesonline.org/fasttree
    # (required by antismash)
    os.system('wget -q https://microbesonline.org/fasttree/FastTreeMP')
    os.system('chmod a+x FastTreeMP')
    os.system('ln -sf FastTreeMP fasttree')
    os.system('ln -sf FastTreeMP FastTree')
    # (tmhmm) tmhmm 2.0c https://services.healthtech.dtu.dk/software.php
    if path_to_tmhmm:
        os.chdir(bin_dir)
        os.system(f'cp {path_to_tmhmm} {bin_dir}')
        os.system(f'tar -xf {path_to_tmhmm.name}')
        perl_exec = which('perl')  # \/usr\/bin\/perl on archlinux
        print(f'detected perl at {perl_exec}')
        os.system(f'sed -i "s|/usr/local/bin/perl|{perl_exec}|g" tmhmm-2.0c/bin/tmhmm')
        os.system(f'sed -i "s|/usr/local/bin/perl|{perl_exec}|g" tmhmm-2.0c/bin/tmhmmformat.pl')
        os.system('chmod a+x tmhmm-2.0c/bin/tmhmm')
        os.system('rm tmhmm-2.0c.Linux.tar.gz')
        os.system('ln -s tmhmm-2.0c tmhmm')
    # (signalp) 6.0g
    if path_to_signalp:
        model_dir = Path(
            which('python')).parent.parent / 'lib' / 'python3.11' / 'site-packages' / 'signalp' / 'model_weights'
        print(f'signalp models will be installed at {model_dir}')
        os.chdir(bin_dir)
        os.system(f'cp {path_to_signalp} .')
        os.system(f'tar -xf {path_to_signalp.name}')
        os.system(f'{Path(which("python")).parent / "pip"} install signalp6_fast/signalp-6-package')
        os.system(f'cp -r signalp6_fast/signalp-6-package/models/* {model_dir}')
        os.system(f'rm -rf {path_to_signalp.name} signalp6_fast')

    #(antismash)
    antismash_database_python_dir = Path(
        which('python')).parent.parent / 'lib' / 'python3.11' / 'site-packages' / 'antismash' / 'databases'
    os.chdir(bin_dir)
    os.system('wget https://dl.secondarymetabolites.org/releases/7.0.0/antismash-7.0.0.tar.gz')
    os.system('tar -xf antismash-7.0.0.tar.gz')
    os.system(f'{Path(which("python")).parent / "pip"} install --upgrade ./antismash-7.0.0')
    os.system('rm -rf antismash-7.0.0.tar.gz antismash-7.0.0')
    os.system(f'rm -rf {antismash_database_python_dir}')
    os.system(f'ln -s {path_to_antismash_db} {antismash_database_python_dir}')
