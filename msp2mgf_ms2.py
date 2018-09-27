import argparse
import re
import ntpath

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--in", required=True,
	help="name of the input file")
ap.add_argument("-o", "--out", required=False,
    help="base name of the output file")
ap.add_argument("-m", "--max", required=False,
    help="max size in MB of the output mgf file before it is split - default is 500")
args = vars(ap.parse_args())

if args['out'] is None:
    raw_filename = ntpath.basename(args['in'])
    base_filename = ''.join(raw_filename.split('.')[0:-1])
    args['out'] = base_filename

if args['max'] is None:
    args['max'] = 500

peptide = ''
charge = ''
mw = ''
comment = ''
rawname = ''
mz = []
intensity = []
mgf_outfile = open(args['out']+'.mgf', 'w')
#ms2_outfile = open('./#ms2/'+args['out']+'.#ms2', 'w')
#db_outfile = open('./db/'+args['out']+'.csv', 'w')
suffix = 0
scannum = 1
shouldWrite = False

#db_outfile.write('"Scan ID","Sequence","Charge","Precursor m/z","Mods"\n')

with open(args['in']) as f:
    for line in f:
        if line.startswith('Name:'):
            raw_name = line.split(" ")[1].rstrip("\r\n")
            peptide = raw_name.split("/")[0]
            raw_charge = raw_name.split("/")[1].rstrip("\r\n")
            charge = raw_charge.split("_")[0]
        if line.startswith('MW:'):
            mw = line.split(" ")[1]
        if line.startswith('Comment:'):
            comment = " ".join(line.split(" ")[1:])
            header = dict(re.findall(r'(\S+)=(".*?"|\S+)', comment))
            if 'Acetyl' in header['Mods'] or 'Propionamide' in header['Mods'] or 'Carbamyl' in header['Mods'] or len(peptide)<8 or len(peptide)>30:
                shouldWrite = False
            else:
                shouldWrite = True
        if line.startswith('Num peaks:'):
            if shouldWrite:
                mgf_outfile.write('BEGIN IONS\n')
                mgf_outfile.write('TITLE=')
                mgf_outfile.write(peptide+'\n')
                #mgf_outfile.write(args['out']+'.'+str(scannum)+'.'+str(scannum)+'.'+charge+'\n')
                mgf_outfile.write('RAWSCANS=')
                mgf_outfile.write(str(scannum)+'\n')
                mgf_outfile.write('SCANS=')
                mgf_outfile.write(str(scannum)+'\n')
                mgf_outfile.write('RTINSECONDS=')
                mgf_outfile.write(str(scannum)+'\n')
                mgf_outfile.write('PEPMASS=')
                mgf_outfile.write(str(header['Parent'])+'\n')
                mgf_outfile.write('CHARGE=')
                mgf_outfile.write(charge+'+\n')

                #ms2_outfile.write('S\t'+str(scannum)+'\t'+str(scannum)+'\t'+str(header['Parent'])+'\n')
                #ms2_outfile.write('Z\t'+charge+'\t'+str(float(mw) - float(charge) + 1)+'\n')

                #db_outfile.write('"'+str(scannum)+'.'+str(scannum)+'.'+charge+'","'+peptide+'","'+charge+'","'+str(header['Parent'])+'","'+header['Mods']+'"\n')

        if (line or 'x')[0].isdigit():
            if shouldWrite:
                peak_info = line.split("\t")[0:2]

                mgf_outfile.write(" ".join(peak_info)+'\n')

                #ms2_outfile.write(" ".join(peak_info)+'\n')

        if line in ['\n', '\r\n']:
            if shouldWrite:
                mgf_outfile.write('END IONS\n\n')

            if (mgf_outfile.tell()/(1024*1024)) > int(args['max']):
                mgf_outfile.close()
                #ms2_outfile.close()
                #db_outfile.close()
                new_filename = args['out'] + '_' + str(suffix)
                mgf_outfile = open('./mgf/'+new_filename+'.mgf', 'w')
                #ms2_outfile = open('./#ms2/'+new_filename+'.#ms2', 'w')
                #db_outfile = open('./db/'+new_filename+'.csv', 'w')
                #db_outfile.write('"Scan ID","Sequence","Charge","Precursor m/z","Mods"\n')
                suffix += 1

            scannum += 1


mgf_outfile.close()
#ms2_outfile.close()
#db_outfile.close()