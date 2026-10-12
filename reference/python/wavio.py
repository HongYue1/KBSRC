import struct, numpy as np
def read(path):
    b=open(path,'rb').read(); assert b[:4]==b'RIFF' and b[8:12]==b'WAVE'
    i=12; fmt=None; data=None
    while i+8<=len(b):
        cid=b[i:i+4]; sz=struct.unpack('<I',b[i+4:i+8])[0]; body=b[i+8:i+8+sz]
        if cid==b'fmt ': fmt=struct.unpack('<HHIIHH',body[:16]); ext=body
        elif cid==b'data': data=body
        i+=8+sz+(sz&1)
    tag,ch,sr,_,ba,bits=fmt
    if tag==0xFFFE: tag=struct.unpack('<H',ext[24:26])[0]
    if tag==3: x=np.frombuffer(data,'<f8' if bits==64 else '<f4').astype(np.float64)
    elif bits==16: x=np.frombuffer(data,'<i2')/32768.0
    elif bits==32: x=np.frombuffer(data,'<i4')/2147483648.0
    elif bits==24:
        a=np.frombuffer(data,np.uint8).reshape(-1,3).astype(np.int32)
        v=a[:,0]|(a[:,1]<<8)|(a[:,2]<<16); v=np.where(v>=1<<23,v-(1<<24),v); x=v/8388608.0
    else: raise ValueError(bits)
    return x.reshape(-1,ch), sr, (tag,bits)
def write(path,x,sr,bits=64):
    x=np.asarray(x,np.float64); x=x.reshape(len(x),-1); ch=x.shape[1]
    d=x.astype('<f8' if bits==64 else '<f4').tobytes(); bps=bits//8
    fmt=struct.pack('<HHIIHHH',3,ch,sr,sr*ch*bps,ch*bps,bits,0)
    fact=struct.pack('<I',len(x))
    out=b'WAVE'+b'fmt '+struct.pack('<I',len(fmt))+fmt+b'fact'+struct.pack('<I',4)+fact+b'data'+struct.pack('<I',len(d))+d
    open(path,'wb').write(b'RIFF'+struct.pack('<I',len(out))+out)

