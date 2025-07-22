function ps = plotPeriodLength(peakMatrix, FramToMin)

p1 = peakMatrix(:,1:size(peakMatrix,2)-1);
p2 = peakMatrix(:,2:size(peakMatrix,2));
ps = (p2-p1)*FramToMin;

figure();
x = [];
y = [];
for c = 1:size(ps,2)
    x = [x c*ones(1,size(peakMatrix,1))];
    y = [y ps(:,c)'];
end
swarmchart(x,y);
xlabel('Cycle Number');
ylabel('Period length (min)');

end