#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
#include<Python.h>

float checking_whether_densed_to_next_layer(float x,float y,bool arg_x,bool arg_y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,int sequence_number_in_next_layer_file,int current_layer)
{
	float A ;

	x = 2 * x - arg_x ;
	y = 2 * y - arg_y ;
	
	int arg_x_plus_double_arg_y = 2 * arg_y + arg_x ;

	int arg = layer_length_accumulated_in_front_of_each_layer[current_layer] + 4 * sequence_number_in_next_layer_file + arg_x_plus_double_arg_y ;
	
	if ( all_layer_whether_densed[arg] )
	{
		arg_x = int( 2 * x ) ;
		arg_y = int( 2 * y ) ;

		sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file[arg] ;

		A = checking_whether_densed_to_next_layer(x,y,arg_x,arg_y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,sequence_number_in_next_layer_file,current_layer+1) ;
	}
	else
	{
		float A1 , A2 , A3 , A4 ;
		A1 = all_layer_corner_mag[ 4 * arg ] ;
		A2 = all_layer_corner_mag[ 4 * arg + 1 ] ;
		A3 = all_layer_corner_mag[ 4 * arg + 2 ] ;
		A4 = all_layer_corner_mag[ 4 * arg + 3 ] ;
		
		A = (1-x)*(1-y)*A1 + x*(1-y)*A2 + (1-x)*y*A3 + x*y*A4 ;
	}

	return A ;
}

float interpolating_to_get_magnification(float x,float y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
	if( x<=-box_size || x>=box_size || y<=-box_size || y>=box_size )
        {
            float s = pow(x*x+y*y,0.5);
            return (s*s+2) / s / pow(s*s + 4,0.5);
        }
 	
	float A ;
	
	x = 0.5 * ( x + box_size ) / box_size ;
	y = 0.5 * ( y + box_size ) / box_size ;

	if ( all_layer_whether_densed[0] )
	{
		bool arg_x = int( 2 * x ) ;
		bool arg_y = int( 2 * y ) ;

		int sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file[ layer_length_accumulated_in_front_of_each_layer[0] + 0 * 4 + 0 ] ;

		A = checking_whether_densed_to_next_layer(x,y,arg_x,arg_y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,sequence_number_in_next_layer_file,0+1) ;
	}
	else
	{
		float A1 , A2 , A3 , A4 ;
                A1 = all_layer_corner_mag[ 0 ] ;
                A2 = all_layer_corner_mag[ 1 ] ;
                A3 = all_layer_corner_mag[ 2 ] ;
                A4 = all_layer_corner_mag[ 3 ] ;

                A = (1-x)*(1-y)*A1 + x*(1-y)*A2 + (1-x)*y*A3 + x*y*A4 ;

	}
	
	return A ;
} 

void generating_magnification_lightcurve(float *hjds,float *lightcurve,float *xtraj_list,float *ytraj_list,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        float t0,u0,te,alpha;
        t0 = parmfit[0];
        u0 = parmfit[1];
        te = parmfit[2];
        alpha = parmfit[3];

        float cos_alpha,sin_alpha;
        cos_alpha = cos(alpha/180*3.1415926);
        sin_alpha = sin(alpha/180*3.1415926);

	float xtraj , ytraj ;

	for(int i=0;i<nhjd;i++)
        {
                xtraj = cos_alpha * (hjds[i]-t0)/te - u0*sin_alpha ;
                ytraj = sin_alpha * (hjds[i]-t0)/te + u0*cos_alpha ;
		
		xtraj_list[i] = xtraj;
		ytraj_list[i] = ytraj;

		lightcurve[i] = interpolating_to_get_magnification(xtraj,ytraj,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
	}

}

float getchi2_sub(float *iflux, float *iferr, float *lc, int nhjd)
{
        double d0=0,d1=0,b00=0,b01=0,b10=0,b11=0;
        double wght,det;
        double fs,fb;
        for(int i=0;i<nhjd;i++)
        {
                wght = iflux[i]/iferr[i]/iferr[i];
                d0 += wght*lc[i];
                d1 += wght;
                b00 += lc[i]*lc[i]/iferr[i]/iferr[i];
                b01 += lc[i]/iferr[i]/iferr[i];
                b11 += 1/iferr[i]/iferr[i];
        }
        b10 = b01;
        //std::cout<<b00<<" "<<b01<<" "<<b11<<std::endl;

        fs = b11*d0-b01*d1;
        fb = -b10*d0+b00*d1;
        det = b00*b11-b10*b01;
        if(det != 0)
        {
                fs /= det;
                fb /= det;
        }
        else
        {
                return -1;
                //fs = 0;
                //fb = d0/b11;
        }
        float res,chi2=0;
        for(int i=0;i<nhjd;i++)
        {
                //std::cout<<(iflux[i]-fb)/fs<<" "<<lc[i]<<std::endl;
                res = iflux[i]-lc[i]*fs-fb;
                //cout << res << endl;
                chi2 += res*res/iferr[i]/iferr[i];
        }
	return chi2;
}

float getchi2(float *hjds,float *iflux, float *iferr,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        float *lightcurve = new float[nhjd];
	float *xtraj_list = new float[nhjd];
	float *ytraj_list = new float[nhjd];
        float chi2;

        generating_magnification_lightcurve(hjds,lightcurve,xtraj_list,ytraj_list,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);

        chi2 = getchi2_sub(iflux,iferr,lightcurve,nhjd);
        
	delete [] lightcurve;
        delete [] xtraj_list;
	delete [] ytraj_list;


	return chi2;
}

int printlc(float *hjds,int nhjd, float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        float *lc = new float[nhjd];
        float *xtraj_list = new float[nhjd];
        float *ytraj_list = new float[nhjd];

        generating_magnification_lightcurve(hjds,lc,xtraj_list,ytraj_list,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
        //std::ofstream outFile("lc.txt");
        //for(int i=0;i<nhjd;i++)
        //{
        //        outFile<<hjds[i]<<" "<<lc[i]<<" "<<xtraj_list[i]<<" "<<ytraj_list[i]<<" "<<parmfit[0]<<" "<<parmfit[1]<<" "<<parmfit[2]<<" "<<parmfit[3]<<" "<<std::endl;
        //}
        //outFile.close();
        delete [] lc;
        delete [] xtraj_list;
        delete [] ytraj_list;

	return 1;

}

extern "C"
{

void * wrapgetchi2(float *hjds,float *iflux, float *iferr,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size,float *chi2)
{
        *chi2 = getchi2(hjds,iflux,iferr,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
}

void * wrapprintlc(float *hjds,int nhjd,float *parmfit,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size)
{
        printlc(hjds,nhjd,parmfit,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
}

void * wrapinterpolating_to_get_magnification(float x,float y,float *all_layer_corner_mag,bool *all_layer_whether_densed,int *all_layer_sequence_number_in_next_layer_file,int *layer_length_accumulated_in_front_of_each_layer,float box_size,float *A)
{
        *A = interpolating_to_get_magnification(x,y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size);
}

}











