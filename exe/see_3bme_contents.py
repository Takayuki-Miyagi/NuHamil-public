#!/usr/bin/env python3
"""
You can use this to check a *.stream.bin file.
Usage:
  python3 see_file_contents.py hoge0 hoge1 hoge2 hoge3
  hoge0: filename
  hoge1: byte width of entries
  hoge2: number of first lines you want to see. ex) 1, 2, 5, 100
         number of last lines you want to see. ex) -1, -2, -5, -100
         first and final entries you want to see. ex) 1-10, 3432-5445
"""
import os, sys
import numpy as np

def see_header(fn, elm_size=4, line=10):
  f = open(fn,'rb')
  for i in range(line):
    if(elm_size==2): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float16)
    if(elm_size==4): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float32)
    if(elm_size==8): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float64)
    print("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}".format(*arr))
  f.close()

def see_tail(fn, elm_size=4, line=10):
  f = open(fn,'rb')
  f.seek(0,os.SEEK_END)
  n_elm = int(f.tell() / elm_size)
  n_rest = n_elm%10
  f.seek(( -10*(line-1)-n_rest) * elm_size,os.SEEK_END)
  for i in range(line-1):
    if(elm_size==2): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float16)
    if(elm_size==4): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float32)
    if(elm_size==8): arr = np.frombuffer(f.read(elm_size*10),dtype=np.float64)
    print("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}".format(*arr))
  if(n_rest > 0):
    if(elm_size==2): arr = np.frombuffer(f.read(elm_size*n_rest),dtype=np.float16)
    if(elm_size==4): arr = np.frombuffer(f.read(elm_size*n_rest),dtype=np.float32)
    if(elm_size==8): arr = np.frombuffer(f.read(elm_size*n_rest),dtype=np.float64)
    print("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}".format(*arr))
  f.close()

def see_middle(fn, elm_size=4, i_begin=0, i_end=100):
  f = open(fn,'rb')
  f.seek((i_begin-1)*elm_size)
  for i in range(i_begin,i_end+1):
    if(elm_size==2): arr = np.frombuffer(f.read(elm_size),dtype=np.float16)
    if(elm_size==4): arr = np.frombuffer(f.read(elm_size),dtype=np.float32)
    if(elm_size==8): arr = np.frombuffer(f.read(elm_size),dtype=np.float64)
    print("{:18d} {:12.6f}".format(i,arr[0]))
  f.close()

if(__name__=="__main__"):
  if(len(sys.argv)==1):
    print("See usage!")
  elif(len(sys.argv)==2):
    see_header(sys.argv[1])
  elif(len(sys.argv)==3):
    see_header(sys.argv[1], int(sys.argv[2]))
    print()
  elif(len(sys.argv)==4):
    string = sys.argv[3]
    strings = sys.argv[3].split("-")
    if( len(strings) == 1 ):
      line = int(string)
      see_header(sys.argv[1], int(sys.argv[2]), line=line)
    elif( len(strings) == 2 ):
      if(strings[0]==""):
        line = int(string)
        see_tail(sys.argv[1], int(sys.argv[2]), line=-line)
      else:
        i_begin = int(strings[0])
        i_end = int(strings[1])
        see_middle(sys.argv[1], int(sys.argv[2]), i_begin=i_begin, i_end=i_end)
    else:
      print("See usage!")
  else:
    print("See usage!")
