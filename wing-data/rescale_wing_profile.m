clear, close all, clc

filename = 'naca0008_CL08_NP628e0.dat';

data = importdata(filename);

figure(1); plot(data(2:end,1),data(2:end,2)); hold on;

xmax = max(data(2:end,1));

data(2:end,:) = data(2:end,:)/xmax;

plot(data(2:end,1),data(2:end,2));

fid =fopen([filename(1:end-4),'_adimensional.dat'],'w');
fprintf(fid,'%d\t',[data(1,1) data(1,2)]);
fprintf(fid,'\n');
for i = 2:length(data(:,1))
  fprintf(fid,'%e\t',data(i,:));
  fprintf(fid,'\n');
end
fclose(fid);
	
