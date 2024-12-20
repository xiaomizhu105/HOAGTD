  

    nSig  = 30;
    inSig = nSig;
    nSig = nSig/255;
    O_Img = double(imread('Monarch.tif'))./255;

    randn('seed', 0);
    N_Img = O_Img + nSig* randn(size(O_Img));                                   %Generate noisy image
    PSNR  =  csnr( N_Img, O_Img, 0, 0 );
    fprintf( 'Noisy Image: nSig = %2.3f, PSNR = %2.2f \n\n\n', inSig, PSNR );
    
    Par   = ParSet(nSig);   
    E_Img = WSNM_DeNoising( N_Img, O_Img, Par );                                %WNNM denoisng function
    PSNR  = csnr( O_Img, E_Img, 0, 0 );
    
    fprintf( 'Estimated Image: nSig = %2.3f, PSNR = %2.2f \n\n\n', nSig, PSNR );
    imshow(uint8(E_Img.*255));