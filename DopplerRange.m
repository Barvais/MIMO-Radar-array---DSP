%Performing the Doppler Range FFT, with opption to print results.
%TR - Post Range FFT Matrix.
%R_max - Maximum range as calculated.
%Print - Prints results if recives 1.
function tr = DopplerRange(TR,R_max,Print)
    
    r1 = linspace(0,R_max,length(TR));
    A = size(TR);
    sp = -floor(A(2)/2):(floor(A(2)/2)-1);
    tr = fftshift(fft(TR,length(sp),2),2);
   
    if(Print == 1)
        figure
        trAbs = 10*log10(abs(tr));
        imagesc(sp,r1,trAbs)
        colorbar
        xlim([-10 10])
        title("Doppler-Range estimation");
        ylabel("Range [m]")
        xlabel("Doppler Velocity [m/s]")
    end
end

