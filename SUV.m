% load data & DICOM infos
PET = read_folder(folder(j).name);
tmp = dir(folder(j).name);
temp = dicominfo([folder(j).name,'\',tmp(3).name]);
clear tmp
% PET DICOM scaling
PET = PET * temp.RescaleSlope + temp.RescaleIntercept;
% SUV calculation
Dose            = temp.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;
InjTime         = temp.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
AcqTime         = temp.AcquisitionTime;
TimeDif         = ( str2num(AcqTime(1:2))*3600 + str2num(AcqTime(3:4))*60 + str2num(AcqTime(5:6)) )  -  ( str2num(InjTime(1:2))*3600 + str2num(InjTime(3:4))*60 + str2num(InjTime(5:6)) );
weight          = temp.PatientWeight;
HalfLife        = temp.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;
actDose         = Dose * 0.5^(TimeDif/HalfLife);
PET    = 1000*PET*weight/actDose;
clear Dose InjTime AcqTime TimeDif weight HalfLife actDose temp
