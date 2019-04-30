#!/usr/bin/python

import csv
import operator
import sys
import numpy as np

# dictionary mapping row[0] to a list of (int(row[1]), row) values
report_map_CG_Plus = {}
report_map_CG_Minus = {}
report_map_CHG_Plus = {}
report_map_CHG_Minus = {}
report_map_CHH_Plus = {}
report_map_CHH_Minus = {}

with open("../BS_Forward.fastq_bismark_pe.CX_report.txt", 'r', newline='') as report:
    reader = csv.reader(report, delimiter='\t' )
    for row in reader:
        if row[2] != "+" or row[5] != "CG":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CG_Plus.setdefault(key, []).append((value, line))
		
		if row[2] != "-" or row[5] != "CG":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CG_Minus.setdefault(key, []).append((value, line))
		
		if row[2] != "+" or row[5] != "CHG":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CHG_Plus.setdefault(key, []).append((value, line))
		
		if row[2] != "-" or row[5] != "CHG":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CHG_Minus.setdefault(key, []).append((value, line))
		
		if row[2] != "+" or row[5] != "CHH":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CHH_Plus.setdefault(key, []).append((value, line))
		
		if row[2] != "-" or row[5] != "CHH":
            continue
        key, value = row[0], int(row[1]) 
        line = '    '.join(row)
        report_map_CHH_Minus.setdefault(key, []).append((value, line))
        
with open("../../Reapeat_analysis/FIMO/BaSAT2/fimo_out/fimo_BaSAT2.gff", 'r', newline='') as fimo, \
     open("x_CG+.txt", 'w', newline='') as fout_CG_Plus, \
	 open("x_CG-.txt", 'w', newline='') as fout_CG_Minus, \
	 open("x_CHG+.txt", 'w', newline='') as fout_CHG_Plus, \
	 open("x_CHG-.txt", 'w', newline='') as fout_CHG_Minus, \
	 open("x_CHH+.txt", 'w', newline='') as fout_CHH_Plus, \
	 open("x_CHH-.txt", 'w', newline='') as fout_CHH_Minus:
	 
    reader = csv.reader(fimo, delimiter='\t')
    writer = csv.writer(fout_CG_Plus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])
		
		
	writer = csv.writer(fout_CG_Minus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])


	writer = csv.writer(fout_CHG_Plus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])


	writer = csv.writer(fout_CHG_Minus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])


	writer = csv.writer(fout_CHH_Plus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])


	writer = csv.writer(fout_CHH_Minus, delimiter="\t")
    for row in reader:
        s0 = row[0]
        s1, s2 = map(int, row[3:5])
        if s0 not in report_map:
            continue
        lt = [r for i, r in report_map[s0] if s1 <= i <= s2]
        if len(lt) <= 0:
            continue
        mt = []
        pt = []
        count = 0
        c = 0
        d = 0
        for a in lt:
            a = a.split(",")
            for b in a:
                b = b.split("   ")
                c = c + int(b[3])
                d = d + int(b[4])
                if int(b[3])>0:
                    count = count + 1
                if (int(b[3])+int(b[4])) >= 4:
                    pt.append(str((float(b[3])/((float(b[3])+float(b[4]))))*100))
                    
                else:
                    pt.append(str("Read less than 4"))
                    
                mt.append(b[3:5])
        try:
            z = (float(c)/(float(c)+float(d)))*100
        except ZeroDivisionError:
            z = 0
        pt_m = 0 
        counter  = 0
        for e in pt:
            if e != "Read less than 4":
                pt_m = pt_m + float(e)
                counter += 1
        try:        
            pt_mean = float(pt_m)/float(counter)
        except ZeroDivisionError:
            pt_mean = 0

        writer.writerow(["{}-{}:{}".format(s1, s2, s0), len(lt), count, (int(count)/int(len(lt))*100), c, d, z,  pt_mean, pt, (mt)])

