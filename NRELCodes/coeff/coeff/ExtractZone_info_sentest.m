function Zoneinfo = ExtractZone_info_sentest(data)

for i =2:length(data)
    if strcmp(data{1,i}(1),{'tout'})
        Tamb=data{1,i}(2:end);
        Zoneinfo.Tamb=str2double(Tamb);
        
    elseif strcmp(data{1,i}(1),'s')
        Sol_rad=data{1,i}(2:end);
        Zoneinfo.Sol_rad=str2double(Sol_rad);
        
    elseif strcmp(data{1,i}(1),'i')
        Qint=data{1,i}(2:end);
        Zoneinfo.Qint=str2double(Qint);
        
    elseif strcmp(data{1,i}(1),'t1')
        t1=data{1,i}(2:end);
        Zoneinfo.t1=str2double(t1);
        
    elseif strcmp(data{1,i}(1),'t2')
        t2=data{1,i}(2:end);
        Zoneinfo.t2=str2double(t2);
        
    elseif strcmp(data{1,i}(1),'t3')
        t3=data{1,i}(2:end);
        Zoneinfo.t3=str2double(t3);
        
    elseif strcmp(data{1,i}(1),'m')
        mf=data{1,i}(2:end);
        Zoneinfo.mf=str2double(mf);
        
    elseif strcmp(data{1,i}(1),'rh')
        rh=data{1,i}(2:end);
        Zoneinfo.rh=str2double(rh);
        
    elseif strcmp(data{1,i}(1),'sp')
        sp=data{1,i}(2:end);
        Zoneinfo.sp=str2double(sp);
        
    elseif strcmp(data{1,i}(1),'sp0')
        sp0=data{1,i}(2:end);
        Zoneinfo.sp0=str2double(sp0);
        
%     elseif strcmp(data{1,i}(1),'m1')
%         m1=data{1,i}(2:end);
%         Zoneinfo.m1=str2double(m1);
%         
%     elseif strcmp(data{1,i}(1),'rh1')
%         rh1=data{1,i}(2:end);
%         Zoneinfo.rh1=str2double(rh1);
%         
%     elseif strcmp(data{1,i}(1),'csp1')
%         csp1=data{1,i}(2:end);
%         Zoneinfo.csp1=str2double(csp1);
%         
%     elseif strcmp(data{1,i}(1),'m2')
%         m2=data{1,i}(2:end);
%         Zoneinfo.m2=str2double(m2);
%         
%     elseif strcmp(data{1,i}(1),'rh2')
%         rh2=data{1,i}(2:end);
%         Zoneinfo.rh2=str2double(rh2);
%         
%     elseif strcmp(data{1,i}(1),'csp2')
%         csp2=data{1,i}(2:end);
%         Zoneinfo.csp2=str2double(csp2);
%         
%     elseif strcmp(data{1,i}(1),'m3')
%         m3=data{1,i}(2:end);
%         Zoneinfo.m3=str2double(m3);
%         
%     elseif strcmp(data{1,i}(1),'rh3')
%         rh3=data{1,i}(2:end);
%         Zoneinfo.rh3=str2double(rh3);
%         
%     elseif strcmp(data{1,i}(1),'csp3')
%         csp3=data{1,i}(2:end);
%         Zoneinfo.csp3=str2double(csp3);
%         
%     elseif strcmp(data{1,i}(1),'m4')
%         m4=data{1,i}(2:end);
%         Zoneinfo.m4=str2double(m4);
%         
%     elseif strcmp(data{1,i}(1),'rh4')
%         rh4=data{1,i}(2:end);
%         Zoneinfo.rh4=str2double(rh4);
%         
%     elseif strcmp(data{1,i}(1),'csp4')
%         csp4=data{1,i}(2:end);
%         Zoneinfo.csp4=str2double(csp4);
        
    end
end