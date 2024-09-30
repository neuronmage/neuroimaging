import os

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return (template, outtype, annotation_classes)

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    
    allowed template fields - follow python string module: 
    
    item: index within category 
    subject: participant id 
    seqitem: run number during scanning
    subindex: sub index within group
    """
    
    rs = create_key('rsfmri/rest_run{item:03d}/rest', outtype=('dicom', 'nii.gz'))
    dwi = create_key('dmri/dwi_{item:03d}', outtype=('dicom', 'nii.gz'))
    t1 = create_key('anatomy/T1_{item:03d}', outtype=('dicom', 'nii.gz'))
    bold = create_key('bold/bold_{item:03d}/bold', outtype=('dicom', 'nii.gz'))
    info = {rs: [], dwi: [], t1: [], bold: []}
    last_run = len(seqinfo)
    for s in seqinfo:
        x,y,sl,nt = (s[6], s[7], s[8], s[9])
        if (sl == 186) and (nt == 1) and ('T1' in s[12]):
            info[t1].append(s[2])
        elif (nt == 180) and ('Resting' in s[12]):
            info[rs].append(int(s[2]))
        elif (nt == 200) and ('MOT' in s[12]):
            info[bold].append(s[2])
            last_run = s[2]
        elif (sl > 1) and ('DTI' in s[12]):
            info[dwi].append(s[2])
        elif ('field_mapping' in s[12]):
            info[fm1].append(s[2])
        else:
            pass
    return info
