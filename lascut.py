import laspy
import sys	
import random
import math
import os

def main():
	tt = (35.7589,139.563077) #musashino
	num = 9 # 直交座標原点系
	wx = 250 # 経度方向長さ(m)
	wy = 210 # 緯度方向長さ
	ratio = 0.4 # 間引き率(1=で間引きなし)
	outf = "musashino.txt" # 出力ファイル
	datapath = "data/" 

	def pt2fn(num,px,py):
		pxh = math.floor(px/100)
		pyh = math.floor(py/100)
		fn = "0"+str(num)+chr(65+pyh)+chr(65+pxh)+ \
			str(math.floor(py/10)-pyh*10)+str(math.floor(px/10)-pxh*10)+str(py%10)+str(px%10) 
		return fn 
	
	tilex = 400 
	tiley = 300 
	l2 = latlng2xy()
	ret = l2.proceed(num,tt[0],tt[1])
	px = math.floor((ret[0]+160000)/tilex) 
	dx = (ret[0]%tilex)
	py = math.floor(-(ret[1]-300000)/tiley)
	dy = (ret[1]%tiley)
	fn = pt2fn(num,px,py)
	files = [datapath+fn+".las"]

	wxr = wx/2
	wyr = wy/2
	sx = dx - wxr
	ex = dx + wxr 
	sy = dy - wyr
	ey = dy + wyr	

	if(sx<0):
		files.append(datapath+pt2fn(num,px-1,py)+".las")
	if(sx<0 and sy<0):
		files.append(datapath+pt2fn(num,px-1,py+1)+".las")
	if(sx<0 and ey>=tiley):
		files.append(datapath+pt2fn(num,px-1,py-1)+".las")
	if(sy<0):
		files.append(datapath+pt2fn(num,px,py+1)+".las")
	if(ey>=tiley):
		files.append(datapath+pt2fn(num,px,py-1)+".las")
	if(ex>=tilex):
		files.append(datapath+pt2fn(num,px+1,py)+".las")
	if(ex>=tilex and sy<0):
		files.append(datapath+pt2fn(num,px+1,py+1)+".las")		
	if(ex>=tilex and ey>=tiley):
		files.append(datapath+pt2fn(num,px+1,py-1)+".las")

	data = {
		"ratio": ratio, 
		"latlng": (tt[0],tt[1]),
		"center": (ret[0],ret[1]),
		"xa": (sx,ex),
		"ya": (sy,ey),
		"files": files,
		"out": outf
	}
	print(data) 

	random.seed() 
	flag = 1 
	for f in data['files']:
		if not os.path.isfile(f):
			print("file not found "+f)
			flag = 0
	if(not flag): return 

	of = open(data['out'],"w")
	l = GetLas()
	for f in data['files']:
		l.proceed(f,of,data)
	of.close()
	print(f"total {GetLas.count} pts") 
	print("complete")
	

class GetLas:
	first = True 
	count = 0 

	def printdata(self,cls,x,y,z,r,g,b,out):
			x = round(x - self.center[0] + 0.005,2)
			y = round(y - self.center[1] + 0.005,2)
			z = round(z - self.zm[0] + 0.005,2)
			r = r//256 
			g = g//256
			b = b//256
			print(f"{cls},{x},{y},{z},{r},{g},{b}",file=out)			
	
	def proceed(self,fn,out,data):
		print("open file "+fn)
		xa = data['xa']
		ya = data['ya']
		ratio = data['ratio']
		las = laspy.read(fn)
		header = las.header

		if GetLas.first:
			self.center = data['center']
			tilez = (math.floor(self.center[0]/400)*400,math.floor(self.center[1]/300)*300)
			self.xmin = xa[0]+tilez[0]
			self.xmax = xa[1]+tilez[0] 
			self.ymin = ya[0]+tilez[1] 
			self.ymax = ya[1]+tilez[1]
			self.xm = (header.x_min,header.x_max)
			self.ym = (header.y_min,header.y_max)
			self.zm = (header.z_min,header.z_max)	
			print(f"{data['latlng'][0]},{data['latlng'][1]},{self.xmin},{self.xmax},{self.ymin},{self.ymax},{self.zm[0]},{self.zm[1]}",file=out)
			GetLas.first = False 
			GetLas.count = 0 
	
		#print(header)
		x_scale = header.x_scale
		x_offset = header.x_offset
		y_scale = header.y_scale
		y_offset = header.y_offset
		z_scale = header.z_scale
		z_offset = header.z_offset
		for pt in las.points.array:
			cls = pt[5]
			xs = pt[0] * x_scale + x_offset
			ys = pt[1] * y_scale + y_offset 
			zs = pt[2] * z_scale + z_offset 

			if xs < self.xmin or xs >=self.xmax or ys < self.ymin or ys >= self.ymax :
				continue 
			if random.random()>ratio: continue 
			self.printdata(cls,xs,ys,zs,pt[10],pt[11],pt[12],out)
			GetLas.count  += 1
			if GetLas.count  % 10000 == 0: print(GetLas.count ) 

# 緯度経度 ->　平面直角座標系 変換
# source https://www.gsi.go.jp/common/000061216.pdf
# 河瀬和重 (2011): Gauss-Krüger投影における経緯度座標及び平面直角座標相互間の座標換算についてのより簡明な計算方法，
class latlng2xy:
	def __init__(self):
		a=6378137 ; rf=298.257222101 ; 
		self.m0=0.9999 ; 
		self.s2r=math.pi/648000 ; 
		self.n=0.5/(rf-0.5)
		n = self.n 
		self.anh=0.5*a/(1+n) ; nsq=n*n
		self.e2n=2*math.sqrt(n)/(1+n) ; 
		self.ra=2*self.anh*self.m0*(1+nsq/4+nsq*nsq/64)
		# 展開パラメータの事前入力
		self.alp=[
			0,
			(1/2+(-2/3+(5/16+(41/180-127/288*n)*n)*n)*n)*n,
			(13/48+(-3/5+(557/1440+281/630*n)*n)*n)*nsq,
			(61/240+(-103/140+15061/26880*n)*n)*n*nsq,
			(49561/161280-179/168*n)*nsq*nsq,
			34729/80640*n*nsq*nsq
		]
		# 平面直角座標の座標系原点の緯度を度単位で、経度を分単位で格納
		self.phi0=[0,33,33,36,33,36,36,36,36,36,40,44,44,44,26,26,26,26,20,26]
		self.lmbd0=[0,7770,7860,7930,8010,8060,8160,8230,8310,8390,8450,8415,8535,8655,8520,7650,7440,7860,8160,9240]
		# 該当緯度の 2 倍角の入力により赤道からの子午線弧長を求める関数
	def Merid(self,phi2):
		jt=5 ; jt2=2*jt ; ep=1.0 ;
		n15=1.5*self.n ; 
		s=[0]*20 ; t=[0]*20 ;	 e=[0]*20
		for k in [1,2,3,4,5]:
			e[k]=n15/k-self.n ;
			ep*=e[k] 
			e[k+jt]=n15/(k+jt)-self.n
		dc=2.0*math.cos(phi2) ; 
		s[1]=math.sin(phi2)
		for i in [1,2,3,4,5,6,7,8,9,10]:
			s[i+1]=dc*s[i]-s[i-1] ; 
			t[i]=(1.0/i-4.0*i)*s[i]
		sum=0.0 ; c1=ep ; j=jt
		while(j>0):
			c2=phi2 ; c3=2.0 ; l=j ; m=0
			while(l>0):
				c3/=e[l]
				c2+=(c3)*t[m+1]
				c3*=e[2*j-(l-1)]
				c2+=(c3)*t[m+2]
				l -= 1 ; m+=2 
			sum+=c1*c1*c2 ; c1/=e[j]
			j -= 1 
		return self.anh*(sum+phi2)

	def proceed(self,num,lng,lat):
		def sinh(x):
			return 0.5*(math.exp(x)-math.exp(-x))
		def cosh(x):
			return 0.5*(math.exp(x)+math.exp(-x))
		def arctanh(x):
			return 0.5*math.log((1+x)/(1-x))

		phirad=lng*3600*self.s2r
		lmbdsec = lat * 3600
		# 実際の計算実行部分
		sphi=math.sin(phirad) ; 
		nphi=(1-self.n)/(1+self.n)*math.tan(phirad)
		dlmbd=(lmbdsec-self.lmbd0[num]*60)*self.s2r
		sdlmbd=math.sin(dlmbd) ; 
		cdlmbd=math.cos(dlmbd)
		tchi=sinh(arctanh(sphi)-self.e2n*arctanh(self.e2n*sphi)) ; 
		cchi=math.sqrt(1+tchi*tchi)
		xi=xip=math.atan2(tchi, cdlmbd) ; 
		eta=etap=arctanh(sdlmbd/cchi) ; 
		for j in [5,4,3,2,1]:
			alsin=self.alp[j]*math.sin(2*j*xip) ; 
			alcos=self.alp[j]*math.cos(2*j*xip)
			xi+=alsin*cosh(2*j*etap) ; 
			eta+=alcos*sinh(2*j*etap)

		x=self.ra*xi-self.m0*self.Merid(2*self.phi0[num]*3600*self.s2r) ; 
		y=self.ra*eta		
		return y,x


main()