
 ALLPATH='D:\Motif_Results\Datasets\ASL\data\';
ALLSubPATH={'God';'I';'Norway';'alive';'all';'answer';'boy';'building';'buy';'change_mind_';'cold';'come';'computer_PC_';'cost';'crazy';'danger';'deaf';'different';'draw';'drink';'eat';'exit';'flashlight';'forget';'girl';'give';'glove';'go';'happy';'head';'hear';'hello';'her';'his_hers';'hot';'how';'hurry';'hurt';'innocent';'is_true_';'joke';'juice';'know';'later';'lose';'love';'make';'man';'maybe';'mine';'money';'more';'name';'no';'notmyproblem';'paper';'pen';'please';'polite';'question';'read';'ready';'research';'responsible';'right';'sad';'same';'science';'share';'shop';'soon';'sorry';'spend';'stubborn';'surprise';'take';'temper';'thank';'think';'tray';'us';'voluntary';'wait_notyet_';'what';'when';'where';'which';'who';'why';'wild';'will';'write';'wrong';'yes';'you';'zero'}


count=1;
for i=1:96%size(ALLPATH,2)
    p=ALLSubPATH{i}
  Path = strcat(  ALLPATH,p,'/');
  Allfiles= dir(Path);
  
  for j=1:size(Allfiles,1)-2
      
      A =csvread([Path,num2str(j),'.csv']);
      
      csvwrite([ALLPATH,num2str(count),'.csv'],A);
      count=count+1;
  end
end