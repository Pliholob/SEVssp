clc
clear
tic
load T_Frame.mat
%% Переменные 

N = 256; % для вычисления 
R = 100;   % доля от 100 интервалов

%% БД для хранения звуков
FileName = {''};
ID_frame = {''};
FileDir = {''};
Frame_Start = [nan];
Frame_End = [nan];
Pd_fr_99 = [nan];
Pd_fr_97 = [nan];
Pd_fr_95 = [nan];
Pd_fr_75 = [nan];
Pd_fr_50 = [nan];
Ex = [nan];
mx = [nan];
NameLet = {''};
TB = table(ID_frame,NameLet,FileName,Frame_Start,Frame_End, Ex, mx, FileDir, ...
Pd_fr_99, Pd_fr_97, Pd_fr_95, Pd_fr_75, Pd_fr_50);
% clear

%% ************************************************************************
 
% disp('*** calculation A0 ***')
% II) Расчет матрицы А для первого частотного интервала
D=fix((N-1)/2);		% количесво границ между частотными интервалами
omega=pi/D;			% ширина первого частотного интервала
% R=(D-1)/2+1;			% количество частотных интервалов
omega = pi/(2*R);			% ширина первого частотного интервала
% disp('*** calculation A ***')
A=zeros(N,N,R);
for i=1:N
for k=1:N
if i==k
A(i,k,1)=(omega)/pi;
else
A(i,k,1)=sin((omega)*(i-k))/(pi*(i-k));
end;
end;
end; 
% На основании матрицы А,
% рассчитываются матрицы для всех других частотных инервалов
for r=2:R
omega_r=2*omega*(r-1); 	% ширина последующих частотных интервалов 
                        % в 2 раза больше ширины 
                        % первого частотного интервала
for i=1:N
for k=1:N
A(i,k,r)=2*A(i,k,1)*cos(omega_r*(i-k));       
end;
end;
end;% for r=2:R
% [Q, L] = eig(As);
%  Qr=fliplr(Qr);       % 
%  Lr=diag(Lr);        %  - выделение диагонали из матрицы   
%  Lr=wrev(Lr);        %  - реверс вектора


% clear i k r D omega_r

%% ------------------------------------------------------------------------- 

% Dir='C:\Users\Kurlova\OneDrive\Project\SSP_SevSU\Later';
currentFolder = pwd;
versn='';
name=['а' 'б' 'в' 'г' 'д' 'е' 'ё' 'ж' 'з' 'и' 'й' 'к' 'л' 'м' 'н' 'о' 'п' 'р' 'с' 'т' 'у' 'ф' 'х' 'ц' 'ч' 'ш' 'щ' 'ы' 'э' 'ю' 'я'];
NBukv=length(name);
dl_bukv=[71 26 17 16 22 19 11 15 21 49 18 29 31 27 33 26 23 30 31 34 18 16 16 20 15 15 12 19 18 11 12];
ext='.wav';
countTF = 0;
for  nbukv=1:NBukv
 for cBukv=1:dl_bukv(nbukv)  
  file=fullfile(currentFolder, '\Later\', name(nbukv), ...
      [name(nbukv) num2str(cBukv) ext versn]);
  [X,Fs]=audioread(file);
  
FileName = [name(nbukv) num2str(cBukv) ext versn];
FileDir = fullfile('\Later\' , name(nbukv), [name(nbukv) num2str(cBukv) ext versn]);

LenFr = length(X);
Frame_Start = 1;

while Frame_Start+N-1 < LenFr
countTF = countTF +1;
Frame_End = Frame_Start+N-1;
ID_frame = [name(nbukv) num2str(cBukv)];
x = X(Frame_Start:Frame_Start+N-1);

mx = mean(x);
Ex = sum(x.^2);
 for r = 1:R
     Pd(r)=(x'*A(:,:,r)*x)/Ex;
 end
 Ps = sort(Pd,'descend');
 Es=0;

P_50 = 0; % !!!!!!!!!!!
P_75 = 0; 
P_99 = 0; 
P_97 = 0; 
P_95 = 0; 
P_h  = 0; 

 for r=2:R
      Es= Es + Ps(r);
if Es <= 0.5 
    P_50 = P_50 +1;
end
if Es <= 0.75 
    P_75 = P_75 +1;
end
if Es <= 0.95 
    P_95 = P_95 +1;
end
if Es <= 0.97 
    P_97 = P_97 +1;
end

if Es <= 0.99 
    P_99 = P_99 +1;
end

% P_h  = 0;
 end

% ------------------------------------------------------------------
    HEx=((2*Df*Ex)/pi);
    r=2;
    mr=r;
    Hm=Ex;    
    for r=2:R
        if Pr(r)<HEx            
       He=abs(Pr(r)-HEx);
        if He<Hm
            mr=r;            
            Hm=He;dolya=Pr(mr)/Ex;
        end;
        end;
    end;
    clear Pr r He Hm HEx


% Заполнение таблицы


 
SR.ID_frame = {[ID_frame '_' num2str(Frame_Start) ...
     '_' num2str(Frame_End)]};
SR.NameLet  = {name(nbukv)};
SR.FileName  = {FileName};
SR.FileDir = {FileDir};
SR.Frame_Start = Frame_Start;
SR.Frame_End = Frame_End;
SR.Ex = Ex;
SR.mx = mx;
SR.Pd_fr_99 = P_99*((Fs/2)/100);
SR.Pd_fr_97 = P_97*((Fs/2)/100);
SR.Pd_fr_95 = P_95*((Fs/2)/100);
SR.Pd_fr_75 = P_75*((Fs/2)/100);
% if P_50 == 0 
%     P_50=1;
% end;
SR.Pd_fr_50 = P_50*((Fs/2)/100);

TB = [TB;struct2table(SR)];
Frame_Start = Frame_Start+N;

 end; % while Frame_Start+N-1 < LenFr
 end; %for cBukv=1:1%dl_bukv(nbukv)  
  file
end; % for  nbukv=1:1%NBukv
TB([1],:) = [];

% save('LaterDS.mat', 'TB') 





















toc