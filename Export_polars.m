wing.polar{1,1}(:,1)=wing.polar{1,1}(:,1)*180/pi;

temp=wing.polar{1,1}(:,3);
wing.polar{1,1}(:,3)=wing.polar{1,1}(:,4);
wing.polar{1,1}(:,4)=temp;

for i=1:stabilizer.n_deflex
    stabilizer.polar{1,i}(:,1)=stabilizer.polar{1,i}(:,1)*180/pi;
    temp=stabilizer.polar{1,i}(:,3);
    stabilizer.polar{1,i}(:,3)=stabilizer.polar{1,i}(:,4);
    stabilizer.polar{1,i}(:,4)=temp;
end

% xlswrite('results1',wing.polar{1,1},'Plan1')
% for i=1:stabilizer.n_deflex
% %   xlswrite('results1',stabilizer.polar{1,i},'Plan2',['A' 7*i-6 ':Z' 7*i])
%     xlswrite('results1',stabilizer.polar{1,i},strcat('Plan',i))
% end    