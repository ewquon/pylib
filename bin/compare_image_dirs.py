#!/usr/bin/env python
"""
Create a beamer presentation comparing images from two different directories
"""
import sys
import os, glob

from beamer import BeamerPresentation

imgtypes = ['png','pdf']

def esc_underscores(s):
    return s.replace('_',r'\_')

if __name__ == '__main__':
    name,dir1,dir2 = sys.argv[1:]
    curdir = os.getcwd()
    dir1 = os.path.join(curdir,dir1)
    dir2 = os.path.join(curdir,dir2)
    dir1name = esc_underscores(os.path.split(dir1)[1])
    dir2name = esc_underscores(os.path.split(dir2)[1])
    fulldir1name = esc_underscores(dir1)
    fulldir2name = esc_underscores(dir2)
    subtitle = fulldir1name + r'\\' + fulldir2name
    imglist = [
        fpath for fpath in glob.glob(os.path.join(dir1,'*'))
        if fpath[-3:] in imgtypes
    ]
    imglist.sort()
    with BeamerPresentation(name,title='Comparison',subtitle=subtitle) as tex:
        for img1path in imglist:
            imgname = os.path.split(img1path)[1]
            img2path = os.path.join(dir2,imgname)
            if not os.path.isfile(img2path):
                print('Skipping',imgname,'because it does not exist in',dir2)
            else:
                frametitle = esc_underscores(os.path.splitext(imgname)[0])
                tex.add_frame_2images(img1path,img2path,
                                      title=frametitle,
                                      subtitles=[dir1name,dir2name])
    os.chdir(name)
    os.system('pdflatex '+name)

