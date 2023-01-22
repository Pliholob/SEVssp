clc
clear
tic
load T_Frame.mat
%% Переменные 

N = 256; % для вычисления 
R = 100;   % доля от 100 интервалов
%% БД для хранения звуков
FileName = {'a.wav';'a.wav'};
ID_frame = [0;0];
Frame_Start = [0;0];
Frame_End = [0;0];
P_99 = [0;0];
P_97 = [0;0];
P_95 = [0;0];
P_h = [0;0];
Ex = [0;0];
mx = [0;0];


%% ************************************************************************
 
disp('*** calculation A0 ***')
% II) Расчет матрицы А для первого частотного интервала
% D=fix((N-1)/2);		% количесво границ между частотными интервалами
% omega=pi/D;			% ширина первого частотного интервала
% R=(D-1)/2+1;			% количество частотных интервалов
omega=2*pi/R;			% ширина первого частотного интервала
disp('*** calculation A ***')
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
clear i k r D omega_r
%%------------------------------------------------------------------------- 

Dir='C:\Users\Kurlova\OneDrive\Project\SSP_SevSU\Later';
versn='';
name=['а' 'б' 'в' 'г' 'д' 'е' 'ё' 'ж' 'з' 'и' 'й' 'к' 'л' 'м' 'н' 'о' 'п' 'р' 'с' 'т' 'у' 'ф' 'х' 'ц' 'ч' 'ш' 'щ' 'ы' 'э' 'ю' 'я'];
NBukv=length(name);
dl_bukv=[71 26 17 16 22 19 11 15 21 49 18 29 31 27 33 26 23 30 31 34 18 16 16 20 15 15 12 19 18 11 12];
ext='.wav';
for  nbukv=1:NBukv
 for cBukv=1:dl_bukv(nbukv)  
  file=fullfile(Dir, name(nbukv), [name(nbukv) num2str(cBukv) ext versn]);
  [X,Fs]=audioread(file);
  
 end;
 file
end;

toc