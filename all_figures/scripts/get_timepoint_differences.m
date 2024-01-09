function DeltaT=get_timepoint_differences
T1="2018-05-29";
T2="2018-10-25";
T3="2018-12-13";
T4="2019-06-04";
T5="2019-10-24";
T6="2022-12-02";

Tall = [T1 T2 T3 T4 T5 T6];
Tall=datetime(Tall);
DeltaT = abs(years(Tall-Tall'));