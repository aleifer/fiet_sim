% simulation script as used in Figure 2 of main paper 

n=50;       %number of neurons
steps=200;  % time-steps in one "iteration". each time-step is ~ 1 burst duration.


beta = .25;	% global inhibition strength
wmax=1;		% single synapse hard bound
m=1;		% Wmax = m*wmax 
pn=2;		 
p=pn/n;		% probability of external stimulation of any neuron at any time
eta=0.125;	% learning rate parameter	
epsilon=0.05;	% strength of heterosynaptic LTD
tau=4;		% time-constant of neural adaptation (only used if alpha is not 0)
alpha =0;	% strength of neural adaptation


%initialize variables
w=zeros(n); dw=zeros(n); dw2=zeros(n); dw3=zeros(n); 
x = zeros(n,1);
y=x;
xdyn=zeros(n,steps);
s=zeros(1,steps);

oldx=x;
oldy=y;
nIters=500; % 10000 in original paper



for iter=1:nIters,
	
	oldw=w;

	for i = 1:steps

		b = rand(n,1)>=(1-p); % set the external stimulation

		%binh = rand(n,1)>=(1-p);            % not relevant
		%y=oldy+1/tau*(-oldy+(oldx)+binh); % y only has to do with neural adaptation.. 
                                          % it can be ignored for cases of
                                          % no adaptation
        y=0;

            % calculate neural activity, x
            x =...                          % neural activity
            ( w*oldx ...                % weights times input frm old activity
            -beta*sum(oldx)...          % global inhibition term
            +b...                       % external stimulation
            -alpha*y)...                % neural adaptation term
            >0;                     % neurons are only active when this 
                                    % sum is greater than zero
                                    
        % STDP 
		dw = eta*(x*double(oldx)'-double(oldx)*x');  %STDP
        
        %Heterosynaptic plasticity
		dw2 = ones(n,1)*max(0, sum(w+dw,1)-m*wmax); %test to see whether sum of weights in each row is too great 
		dw3 = max(0, sum(w+dw,2)-m*wmax)*ones(1,n); %test to see weather sum of weights in each column is too great
        
        %update the new weight matrix
		w=min(wmax,...                  % any weight that exceeds wmax gets set to wmax
            max(0,...                   % force any negative synapses to be zero
            w+dw ...                    % update weight matrix based on STDP
            -eta*epsilon*(dw2+dw3)...   % Heterosynaptic plasticity: penalize all those rows or colums whose sums exceeded the maximum synaptic weight
            -eye(n)*10000*wmax));       % Zero out any synapses from one neuron to itself 
                                        %   (remember that anything
                                        %   negative is set to zero)

		oldx = x;  % remember the old neural activity
		oldy= y;   % remember the old adatation state (not important when no adaption in model)
		xdyn(:,i)=x; % record a copy of the neural activity

        
        %generate a sorted matrix for display purposes
        [~,maxColByRow]=max(w,[],2); %for each row, find the column that has the highest wait
        [~,srtIndx]=sort(maxColByRow); %sort those and record the indices
        
        
    end

        subplot(2,4,1); 
        imagesc(w,[0,wmax]); colormap(hot); colorbar
        title('W'); xlabel('neuron index'); ylabel('neuron index');

        subplot(2,4,2); 
        imagesc(w'*w,[0,wmax^2]); colormap(hot); colorbar
        title('W^T*W'); xlabel('neuron index'); ylabel('neuron index');

        subplot(2,4,3);
        imagesc(w(srtIndx,:))
        title(['w sorted  iter=' num2str(iter)])
        ylabel('DIFFERENT (sorted) neuron index')
        xlabel('neuron index')

        subplot(2,4,4);
        hist(reshape(w,1,[]));
        title('W')
        ylabel('Counts')
        xlabel('Weight')

        
        
        subplot(2,4,5); 
        imagesc(w-oldw); colormap(hot); colorbar
        title('change in W'); xlabel('neuron index'); ylabel('neuron index');


        subplot(2,4,6); 
        imagesc(xdyn); 
        title('neural activity')
        xlabel('time (steps)'); ylabel('neuron index');


        subplot(2,4,7); 
        imagesc(xdyn(srtIndx,:)); 
        title('neural activity')
        xlabel('time (steps)'); ylabel('SORTED neuron index (defined by weights)');



        drawnow; 
    

end






%-----------------------------------------------------
% playback of learned activiy sequences: 
%-----------------------------------------------------
figure;
for iter=1:80,
	
	oldw=w; % not relevant

	for i = 1:steps

		% 4 random bursts to initiate activity: 
		if (iter<5),
		b = rand(n,1)>=(1-p);
		elseif (iter>20)&(iter<25),
		b = rand(n,1)>=(1-p);
		elseif (iter>40)&(iter<45),
		b = rand(n,1)>=(1-p);
		elseif (iter>60)&(iter<65),
		b = rand(n,1)>=(1-p);
		else
		b=0*b; 
        end
        
        
        %binh = rand(n,1)>=(1-p); % not relevant (adaptation)
		%y=oldy+1/tau*(-oldy+(oldx)+binh); % not relevant (adaptation)
        y=0;

		x = (w*oldx-beta*sum(oldx)+b-alpha*y)>0;
		
		oldx = x;
		oldy= y;  
		xdyn(:,i)=x;

	end

	subplot(2,2,1); 
	title(['iter=',num2str(iter)])
	imagesc(w,[0,wmax]); colormap(hot); colorbar
	title('W'); xlabel('neuron index'); ylabel('neuron index');
	subplot(2,2,2); 
	imagesc(w'*w,[0,wmax^2]); colormap(hot); colorbar
	title('W^T*W'); xlabel('neuron index'); ylabel('neuron index');
	subplot(2,2,3); 
	imagesc(xdyn'*xdyn); 
	title(['neural activity, iter =', num2str(iter)]); xlabel('neuron index'); ylabel('neuron index');
	subplot(2,2,4); 
	imagesc(xdyn); 
	title('neural activity')
	xlabel(['time (steps), iter =', num2str(iter)]); ylabel('neuron index');
	drawnow; 
end


