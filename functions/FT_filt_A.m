function [FT_data_filt_A] = FT_filt_A(data, sampFreq)

[b1,a1] = butter(3,1/(sampFreq/2),'high');
f1_data = filtfilt(b1,a1,data')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(3,[119,121]/(sampFreq/2),'stop');
f3_data = filtfilt(b3,a3,f2_data')';

[b4,a4] = butter(8,150/(sampFreq/2),'low');
f4_data = filtfilt(b4,a4,f3_data')';

FT_data_filt_A = f4_data;

end

