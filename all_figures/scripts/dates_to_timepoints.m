function TP = dates_to_timepoints(dates)

%% get sampling TP from dates

dates=string(dates,"yyyy-MM-dd");

TP = zeros(numel(dates),1);
TP(dates=="2018-05-29"|dates=="2018-05-30")=1;
TP(dates=="2018-10-25")=2;
TP(dates=="2018-12-12"|dates=="2018-12-13")=3;
TP(dates=="2019-06-04")=4;
TP(dates=="2019-10-24")=5;
TP(dates=="2022-12-02")=6;

if any(TP==0)
    warning('timepoint not found for at least one date')
end
