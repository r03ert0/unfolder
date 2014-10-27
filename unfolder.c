/*
To compile use:
    gcc -Wall unfolder.c -o unfolder
To use, for example:
./unfolder fibers_thin.trk -braingl -rotate 90 0 -90 -length2width -minMaxLength -cylindre -minMaxLength -filterLength 10 300 -minMaxLength -saveMesh cyl.txt
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Structures
typedef struct
{
	char			id_string[6];
	short			dim[3];
	float			voxel_size[3];
	float			origin[3];
	short			n_scalars;
	char			scalar_name[10][20];
	short			n_properties;
	char			property_name[10][20];
	float			vox_to_ras[4][4];
	char			reserved[444];
	char			voxel_order[4];
	char			pad2[4];
	float			image_orientation_patient[6];
	char			pad1[2];
	unsigned char	invert_x;
	unsigned char	invert_y;
	unsigned char	invert_z;
	unsigned char	invert_xy;
	unsigned char	invert_yz;
	unsigned char	invert_zx;
	int				n_count;
	int				version;
	int				hdr_size;
}TrackHeader;

typedef struct
{
	float	x,y,z;
}float3D;
typedef struct
{
	int	a,b,c;
}int3D;

// Global variables
TrackHeader	hdr;
int			*m;
long		*parr;
float       *sz,SZ;

// Vector algebra
float norm3D(float3D a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
float3D add3D(float3D a, float3D b)
{
	return (float3D){a.x+b.x,a.y+b.y,a.z+b.z};
}
float3D sub3D(float3D a, float3D b)
{
	return (float3D){a.x-b.x,a.y-b.y,a.z-b.z};
}
float3D sca3D(float3D a, float t)
{
	return (float3D){a.x*t,a.y*t,a.z*t};
}
float3D cross3D(float3D a, float3D b)
{
	return (float3D){a.y*b.z-a.z*b.y,b.x*a.z-b.z*a.x,a.x*b.y-a.y*b.x};
}

// Tractography functions
float lengthOfFibre(int j)
{
	float3D	*p;
	int		i;
	float	length;

	length=0;
	p=(float3D*)parr[j];
	for(i=1;i<m[j];i++)
		length+=norm3D(sub3D(p[i],p[i-1]));
	return length;
}
void computeLength(void)		// Compute length
{
	int		j;
	
	for(j=0;j<hdr.n_count;j++)
		printf("%f\n",lengthOfFibre(j));
}
void minMaxLength(void)		// Compute length
{
	int		j;
	float	min,max;
	
	min=max=lengthOfFibre(0);
	for(j=1;j<hdr.n_count;j++)
	{
		if(lengthOfFibre(j)<min)
			min=lengthOfFibre(j);
		if(lengthOfFibre(j)>max)
			max=lengthOfFibre(j);
	}
	printf("length min,max: %f, %f\n",min,max);
}
void lengthHistogram(void)		// Compute length histogram
{
	int		i,j;
	float	length;
	int		nbins=100,*hist;
	float	min=0,max=200;
	
	hist=(int*)calloc(nbins,sizeof(int));
	for(j=0;j<hdr.n_count;j++)
	{
		length=lengthOfFibre(j);		
		i=(int)((nbins-1)*(length-min)/(max-min));
		
		if(i>(nbins-1))
			printf("length: %g\n",length);
	
		hist[i]++;
	}
	for(i=0;i<nbins;i++)
	{
		if(i<nbins-1)
			printf("%i ",hist[i]);
		else
			printf("%i\n",hist[i]);
	}
	free(hist);
}
float3D normal(float3D *p, int np, int i)
{
	float3D	n;
	
	if(i<1)
		n=normal(p,np,i+1);
	else
	if(i==np-1)
		n=normal(p,np,i-1);
	else
	if(i<np-1)
	{
		n=cross3D(sub3D(p[i-1],p[i]),sub3D(p[i+1],p[i]));
		if(norm3D(n)<0.001)
			n=(float3D){0,0,1};
		else
			n=sca3D(n,1/norm3D(n));
	}
	
	return n;
}
void saveMesh(char *path)		// Save as mesh
{
	FILE	*f;
	int		i,j,np,nt;
	float3D	*p,p1,n;
	
	f=fopen(path,"w");
	
	// save header
	np=nt=0;
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		np+=m[j]*2;
		nt+=2*(m[j]-1);
	}
	fprintf(f,"%i %i\n",np,nt);

	if(1)
	{
		// Length-dependent ribbon width
		for(j=0;j<hdr.n_count;j++)
		if(m[j])
		{
			p=(float3D*)parr[j];
			n=(float3D){0,0,0};
			for(i=0;i<m[j];i++)
				n=add3D(n,normal(p,m[j],i));
			n=sca3D(n,sz[j]/norm3D(n));
			for(i=0;i<m[j];i++)
			{
				fprintf(f,"%g %g %g\n",p[i].x,p[i].y,p[i].z);
				p1=add3D(p[i],n);
				fprintf(f,"%g %g %g\n",p1.x,p1.y,p1.z);
			}
		}		
	}
	else
	{
		// Constant ribbon width
		for(j=0;j<hdr.n_count;j++)
		if(m[j])
		{
			p=(float3D*)parr[j];
			n=(float3D){0,0,0};
			for(i=0;i<m[j];i++)
				n=add3D(n,normal(p,m[j],i));
			n=sca3D(n,SZ/norm3D(n));
			for(i=0;i<m[j];i++)
			{
				fprintf(f,"%g %g %g\n",p[i].x,p[i].y,p[i].z);
				p1=add3D(p[i],n);
				fprintf(f,"%g %g %g\n",p1.x,p1.y,p1.z);
			}
		}		
	}
		
	// save triangles
	np=0;
	nt=0;
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];
		for(i=1;i<m[j];i++)
		{
			fprintf(f,"%i %i %i\n",np+2*(i-1),np+2*i,np+2*(i-1)+1);
			fprintf(f,"%i %i %i\n",np+2*(i-1)+1,np+2*i,np+2*i+1);
			nt+=2;
		}
		np+=2*m[j];
	}
	fclose(f);	
}
void mapLengthToSize(void)
{
	int		j;
	float	min,max;
	float	sz0,szmin,szmax;

	szmin=0.1;
	szmax=6;
	min=max=lengthOfFibre(0);
	for(j=1;j<hdr.n_count;j++)
	{
		if(lengthOfFibre(j)<min)
			min=lengthOfFibre(j);
		if(lengthOfFibre(j)>max)
			max=lengthOfFibre(j);
	}

	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		sz0=(lengthOfFibre(j)-min)/(max-min);
		//sz0=pow(sz0,3);
		sz[j]=szmin*(1-sz0)+szmax*sz0;
	}
}
void rotate(float x, float y, float z)
{
    int i,j,n=0;
    float   M[9];
    float3D	*p,p1,p2,origin={0,0,0};
 
    M[0]=cos(z)*cos(y);
    M[1]=-sin(z)*cos(x)+cos(z)*sin(y)*sin(x);
    M[2]=sin(z)*sin(x)+cos(z)*sin(y)*cos(x);
    
    M[3]=sin(z)*cos(y);
    M[4]=cos(z)*cos(x)+sin(z)*sin(y)*sin(x);
    M[5]=-cos(z)*sin(x)+sin(z)*sin(y)*cos(x);
    
    M[6]=-sin(y);
    M[7]=cos(y)*sin(x);
    M[8]=cos(y)*cos(x);
    
	// compute barycentre
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];
		for(i=0;i<m[j];i++)
		{
			origin=add3D(origin,p[i]);
			n++;
		}
	}
	origin=sca3D(origin,1/(float)n);
		
	// apply rotation to vertices
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];

		for(i=0;i<m[j];i++)
		{
			p1=sub3D(p[i],origin);
			p2.x=M[0]*p1.x+M[1]*p1.y+M[2]*p1.z;
			p2.y=M[3]*p1.x+M[4]*p1.y+M[5]*p1.z;
			p2.z=M[6]*p1.x+M[7]*p1.y+M[8]*p1.z;
			p[i]=add3D(p2,origin);
		}
	}		
}
void cylindre(void)
{
	int		i,j,n=0;
	float3D	*p,origin={0,0,0};
	float	r,a,R=0,H=0;
	
	// compute barycentre
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];
		for(i=0;i<m[j];i++)
		{
			origin=add3D(origin,p[i]);
			n++;
		}
	}
	origin=sca3D(origin,1/(float)n);
	
	// compute radius and height
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];
		for(i=0;i<m[j];i++)
		{
			r=sqrt(pow(p[i].x-origin.x,2)+pow(p[i].y-origin.y,2));
			if(r>R)
				R=r;
			if(p[i].z-origin.z>H)
				H=p[i].z-origin.z;
		}
	}
	H=H*2;
	
	// apply cylindrical transformation
	for(j=0;j<hdr.n_count;j++)
	if(m[j])
	{
		p=(float3D*)parr[j];
		for(i=0;i<m[j];i++)
		{
			r=sqrt(pow(p[i].x-origin.x,2)+pow(p[i].y-origin.y,2));
			a=atan2(p[i].y-origin.y,p[i].x-origin.x);
			
			p[i]=(float3D){a*R,r,p[i].z};
		}
	}	
}
int swapint(int i)
{
	unsigned char	*b2;
	b2=(unsigned char *)&i;
	return b2[0]*256*256*256+b2[1]*256*256+b2[2]*256+b2[3];
}
float swapfloat(float f)
{
	unsigned char	b[4],*b2;
	b2=(unsigned char *)&f;
	b[0]=b2[3];
	b[1]=b2[2];
	b[2]=b2[1];
	b[3]=b2[0];
	return *(float*)b;
}

// Main
int main(int argc, char *argv[])
{
	FILE		*f;
	int			i,j,k,n;
	float		*tmp,*pr;
	float3D		*p;
	int			braingl=0;
	
	if(argc<=1)
	{
		printf("Unfolder: a tool to unfold neuroimaging data\n");
		printf("Katja Heuer and Roberto Toro\n");
		printf("Version 1. October 2014\n");
		printf("\n");
		printf("Usage:\n");
		printf("    unfolder data_file [switches...]\n");
		printf("\n");
		printf("Switches are applied and processed in the order (it is not the same to do\n");
		printf("-one -two or -two -one. Available switches:\n");
		printf("    -braingl                             Account for endianness\n");
		printf("                                         idiosyncrasies in brainGL's\n");
		printf("                                         fib to trk converter\n");
		printf("    -rotate x y z                        Rotate the object by x, y and\n");
		printf("                                         then z degrees before unfolding\n");
		printf("    -length2width                        Map fibre length to fibre ribbon\n");
		printf("                                         width\n");
		printf("    -length	                             Display the length of the\n");
		printf("                                         longest fibre\n");
		printf("    -filterVertices min_vertices         Remove all fibres with less than\n");
		printf("                                         min_vertices vertices\n");
		printf("    -width the_width                     Fix the fibre ribbon width to\n");
		printf("                                         the_width\n");
		printf("    -hist                                Display a histogram of fibre\n");
		printf("                                         lengths\n");
		printf("    -minMaxLength                        Compute min and max fibre lengths\n");
		printf("    -cylindre                            Unfold objects using a\n");
		printf("                                         cylindrical transformation\n");
		printf("    -filterLength min_length max_length  Remove all fibres shorter than\n");
		printf("                                         min_length and longer than\n");
		printf("                                         max_length\n");
		printf("    -saveMesh filename.txt               Save the resulting mesh into\n");
		printf("                                         filename.txt\n");
		printf("\n");
		return 0;
	}
// default fiber width
	SZ=0.3;
	
	f=fopen(argv[1],"r");
	
	if(strcmp(argv[2],"-braingl")==0)
		braingl=1;

	// read header
	fread(&hdr,1,sizeof(TrackHeader),f);
	
	if(braingl)
		hdr.n_count=swapint(hdr.n_count);

	m=(int*)calloc(hdr.n_count,sizeof(int));
	sz=(float*)calloc(hdr.n_count,sizeof(float));
	parr=(long*)calloc(hdr.n_count,sizeof(long));
	
	// read tracks
	printf("n_count=%i\n",hdr.n_count);
	pr=(float*)calloc(hdr.n_properties,sizeof(float));
	for(j=0;j<hdr.n_count;j++)
	{
		fread(&n,1,sizeof(int),f);
		if(braingl)
			n=swapint(n);
		tmp=(float*)calloc(n*(3+hdr.n_scalars),sizeof(float));
		fread(tmp,n,(3+hdr.n_scalars)*sizeof(float),f);
		if(braingl)
			for(k=0;k<n*(3+hdr.n_scalars);k++)
				tmp[k]=swapfloat(tmp[k]);
		fread(pr,hdr.n_properties,sizeof(float),f);
		if(braingl)
			for(k=0;k<hdr.n_properties;k++)
				pr[k]=swapfloat(pr[k]);


		m[j]=n;
		p=(float3D*)calloc(m[j],sizeof(float3D));
		parr[j]=(long)p;
		for(i=0;i<m[j];i++)
			p[i]=*(float3D*)&(tmp[i*(3+hdr.n_scalars)]);

		free(tmp);		
	}
	fclose(f);
	free(pr);

	for(i=2+braingl;i<argc;i++)
	{
		printf("%s\n",argv[i]);
		
		if(strcmp(argv[i],"-length")==0)
			computeLength();
		else
		if(strcmp(argv[i],"-width")==0)
		{
			SZ=atof(argv[++i]);
			printf("width: %f\n",SZ);
		}
		else
		if(strcmp(argv[i],"-minMaxLength")==0)
			minMaxLength();
		else
		if(strcmp(argv[i],"-rotate")==0)
		{
			float	x,y,z;
			float	pi=2*acos(0);
			x=atof(argv[++i])*pi/180.0;
			y=atof(argv[++i])*pi/180.0;
			z=atof(argv[++i])*pi/180.0;
			rotate(x,y,z);
		}
		else
		if(strcmp(argv[i],"-length2width")==0)
			mapLengthToSize();
		else
		if(strcmp(argv[i],"-cylindre")==0)
			cylindre();
		else
		if(strcmp(argv[i],"-hist")==0)
			lengthHistogram();
		else
		if(strcmp(argv[i],"-filterVertices")==0)
		{
			int	nmin=atoi(argv[++i]);
			n=0;
			for(j=0;j<hdr.n_count;j++)
				if(m[j]<nmin)
				{
					m[j]=0;
					n++;
				}
			printf("%i filtered out\n",n);
		}
		else
		if(strcmp(argv[i],"-filterLength")==0)
		{
			float lmin=atof(argv[++i]);
			float lmax=atof(argv[++i]);
			n=0;
			for(j=0;j<hdr.n_count;j++)
				if(lengthOfFibre(j)<lmin || lengthOfFibre(j)>lmax)
				{
					m[j]=0;
					n++;
				}
			printf("%i filtered out\n",n);
		}
		else
		if(strcmp(argv[i],"-saveMesh")==0)
			saveMesh(argv[++i]);
	}
	
	return 0;
}