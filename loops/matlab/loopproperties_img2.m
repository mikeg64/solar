%loop coordinates for image sets
npos=12; %number of positions across each coronal loop
nloops=13; %number of loops analysed

tstr171=cell(nloops);
tstr193=cell(nloops);

tsec171=zeros(nloops,1);
tsec193=zeros(nloops,1);

lx171=zeros(nloops,npos);
ly171=zeros(nloops,npos);
lx193=zeros(nloops,npos);
ly193=zeros(nloops,npos);

lxbg193=zeros(nloops,2);
lybg193=zeros(nloops,2);
lxbg171=zeros(nloops,2);
lybg171=zeros(nloops,2);

lheight193=zeros(nloops);
lheight171=zeros(nloops);
llength193=zeros(nloops);
llength171=zeros(nloops);
lwidth193=zeros(nloops);
lwidth171=zeros(nloops);
lddomnode193=zeros(nloops);
lddomnode171=zeros(nloops);
lintensity193=zeros(nloops);
lintensity171=zeros(nloops);

%two points forming a line to sample background intensity
lxbg193(1,:)=[];
lybg193(1,:)=[];
lxbg171(1,:)=[];
lybg171(1,:)=[];
lx171(1,:)=[413 423 435 452 483 515 530 541 557 559 552 529];
ly171(1,:)=[275 230 203 182 164 159 163 171 195 228 267 312];
lx193(1,:)=[420 426 433 442 455 489 509 528 548 558 545 538];
ly193(1,:)=[305 271 238 221 207 183 177 182 200 228 288 312];

%two points forming a line to sample background intensity
 lxbg193(2,:)=[535 547];
 lybg193(2,:)=[203 222];
 lxbg171(2,:)=[546 548];
 lybg171(2,:)=[207 227];
lx171(2,:)=[418 423 433 442 463 490 515 534 557 557 553 534];
ly171(2,:)=[282 235 212 198 179 161 162 168 202 231 261 303];
lx193(2,:)=[430 428 457 480 487 511 539 557 571 575 564 546];
ly193(2,:)=[306 256 212 202 197 191 195 206 227 252 275 307];


%two points forming a line to sample background intensity
 lxbg193(3,:)=[491 500];
 lybg193(3,:)=[226 230];
 lxbg171(3,:)=[543 547];
 lybg171(3,:)=[203 217];
 lx171(3,:)=[416 424 446 455 475 504 530 543 552 554 550 534];
 ly171(3,:)=[273 220 193 188 173 164 176 189 208 246 267 302];
 lx193(3,:)=[432 423 437 451 469 487 499 521 530 537 536 529];
 ly193(3,:)=[309 283 254 256 223 217 219 225 235 251 266 312];
%

%two points forming a line to sample background intensity
 lxbg193(4,:)=[487 497];
 lybg193(4,:)=[228 232];
 lxbg171(4,:)=[515 525];
 lybg171(4,:)=[183 198];
 lx171(4,:)=[416 416 422 426 441 449 461 488 531 546 555 537];
 ly171(4,:)=[297 280 231 217 194 187 182 165 174 195 240 306];
 lx193(4,:)=[430 423 418 421 438 475 500 520 526 530 536 536];
 ly193(4,:)=[316 309 286 256 234 217 217 222 228 234 249 278];
%

%two points forming a line to sample background intensity
% lxbg193(5,:)=[];
% lybg193(5,:)=[];
% lxbg171(5,:)=[];
% lybg171(5,:)=[];
 lx171(5,:)=[413 416 424 441 470 499 522 546 560 557 551 536];
 ly171(5,:)=[267 241 221 200 175 165 168 177 199 227 266 298];
 lx193(5,:)=[424 418 421 430 440 455 478 501 510 520 535 534 ];
 ly193(5,:)=[309 284 269 240 224 215 206 212 215 222 250 280];
%

%two points forming a line to sample background intensity
% lxbg193(6,:)=[];
% lybg193(6,:)=[];
% lxbg171(6,:)=[];
% lybg171(6,:)=[];
 lx171(6,:)=[413 415 448 458 471 509 538 548 557 559 557 531];
 ly171(6,:)=[287 247 194 179 170 150 154 164 184 205 228 313];
 lx193(6,:)=[424 423 427 436 437 447 458 484 512 534 534 521];
 ly193(6,:)=[302 284 254 242 241 231 223 220 228 248 286 321];
%

%after flare

 lx171(7,:)=[425 423 441 449 465 477 495 514 525 536 540 529];
 ly171(7,:)=[305 264 233 219 201 198 194 205 215 231 262 311];
 lx193(7,:)=[427 424 423 429 438 459 474 490 508 532 536 529];
 ly193(7,:)=[317 292 272 241 223 207 206 209 217 234 287 329];

 
 lx171(8,:)=[424 423 436 461 488 515 531 541 541 538 531 523 ];
 ly171(8,:)=[314 274 236 203 193 199 212 235 257 275 296 318];
 lx193(8,:)=[418 410 417 438 454 474 494 514 527 536 539 523];
 ly193(8,:)=[306 278 251 224 212 207 210 219 231 250 274 309];

 lx171(9,:)=[417 421 428 444 477 498 516 532 541 540 537 527];
 ly171(9,:)=[304 255 234 217 195 191 195 213 235 259 284 306];
 lx193(9,:)=[422 419 421 428 437 453 485 499 532 542 540 531];
 ly193(9,:)=[312 290 269 244 227 212 204 202 219 243 269 308];

 lx171(10,:)=[425 423 433 455 474 495 515 530 537 540 537 526];
 ly171(10,:)=[298 252 229 205 194 190 195 208 226 248 274 296];
 lx193(10,:)=[424 420 424 424 429 451 482 500 522 539 544 530];
 ly193(10,:)=[310 295 263 231 217 202 191 193 203 221 254 306];



 lx171(11,:)=[412 416 426 435 446 474 506 520 536 541 538 522];
 ly171(11,:)=[299 271 243 228 216 196 192 201 212 241 270 315];
 lx193(11,:)=[427 427 427 434 438 451 474 501 525 540 547 558 ];
 ly193(11,:)=[314 289 258 228 220 205 197 193 203 224 268 308 ];

lx171(12,:)=[420 427 446 486 523 542 546 530 0 0 0 0 ];
ly171(12,:)=[287 244 213 189 194 218 250 305 0 0 0 0 ];
lx193(12,:)=[413 406 430 466 511 543 554 554 0 0 0 0 ];
ly193(12,:)=[309 247 196 172 166 187 228 292 0 0 0 0 ];

lx171(13,:)=[422 449 467 508 543 562 569 555 0 0 0 0 ];
ly171(13,:)=[260 247 229 212 212 225 255 294 0 0 0 0 ];
lx193(13,:)=[424 450 480 512 536 556 557 547 0 0 0 0 ];
ly193(13,:)=[307 275 251 243 250 261 290 310 0 0 0 0 ];

%lddomnode is the distance of the peak of the loop from the dominant 
%most stable point in the loop system.
%this is the bright point at the node of all the loops at the left hand
%side of the feature under study

tstr171{1}='12_211423';
tsec171(1)=74676;

tstr193{1}='12_211019';
tsec193(1)=73812;

lheight193(1)=152.35;
lheight171(1)=137.14;
llength193(1)=116.9;
llength171(1)=123.11;
lwidth193(1)=3.7;
lwidth171(1)=3.7;
% lintensity193(1)=0;
% lintensity171(1)=0;
 lddomnode193(1)=167.07;
 lddomnode171(1)=179.39;

tstr171{2}='12_214359';
tsec171(2)=76452;
tstr193{2}='12_214043';
tsec193(2)=75636;
 lheight193(2)=130.25;
 lheight171(2)=108.22;
 llength193(2)=108.34;
 llength171(2)=114.84;
 lwidth193(2)=4.9;
 lwidth171(2)=6.1;
%  lintensity193(2)=0;
% lintensity171(2)=0;
  lddomnode193(2)=162.29;
  lddomnode171(2)=173.28;
 
  tstr193{3}='12_221043';
tsec193(3)=77436;
 tstr171{3}='12_221423';
tsec171(3)=78276;
   lheight193(3)=88.68;
  lheight171(3)=127.49;
  llength193(3)=95.05;
  llength171(3)=120.17;
  lwidth193(3)=4.64;
  lwidth171(3)=3.83;
% lintensity193(3)=0;
% lintensity171(3)=0;
  lddomnode193(3)=115.95;
  lddomnode171(3)=170.84;

 tstr193{4}='12_222519';
tsec193(4)=78312;
 tstr171{4}='12_222923';
tsec171(4)=79176;
   lheight193(4)=93.19;
  lheight171(4)=138.01;
  llength193(4)=102.02;
  llength171(4)=124.44;
  lwidth193(4)=5.43;
  lwidth171(4)=6.4;
% lintensity193(4)=0;
% lintensity171(4)=0;
  lddomnode193(4)=123.57;
  lddomnode171(4)=184.63;


tstr193{5}='12_224043';
tsec193(5)=79236;
tstr171{5}='12_224423';
tsec171(5)=80076;
   lheight193(5)=94.02;
  lheight171(5)=120.74;
  llength193(5)=116.37;
  llength171(5)=127.34;
  lwidth193(5)=3.8;
  lwidth171(5)=3.2;
% lintensity193(5)=0;
% lintensity171(5)=0;
  lddomnode193(5)=128.88;
  lddomnode171(5)=158.64;

 tstr193{6}='12_225545';
tsec193(6)=80138;
 tstr171{6}='12_225923';
tsec171(6)=80976;
   lheight193(6)=100.62;
  lheight171(6)=156.6;
  llength193(6)=100.78;
  llength171(6)=115.2;
  lwidth193(6)=4.0;
  lwidth171(6)=4.9;
% lintensity193(6)=0;
% lintensity171(6)=0;
  lddomnode193(6)=121.44;
  lddomnode171(6)=194.74;

%after the flare

 tstr193{7}='12_231008';
tsec193(7)=81001;
 tstr171{7}='12_231423';
tsec171(7)=81876;
   lheight193(7)=110.57;
  lheight171(7)=124.28;
  llength193(7)=102.14;
  llength171(7)=97.72;
  lwidth193(7)=3.6;
  lwidth171(7)=4.8;
% lintensity193(7)=0;
% lintensity171(7)=0;
  lddomnode193(7)=139.78;
  lddomnode171(7)=153.42;

tstr193{8}='12_232443';
tsec193(8)=81876;
tstr171{8}='12_232923';
tsec171(8)=82776;
   lheight193(8)=108.07;
  lheight171(8)=130.81;
  llength193(8)=98.18;
  llength171(8)=91.61;
  lwidth193(8)=6.3;
  lwidth171(8)=6.8;
% lintensity193(8)=0;
% lintensity171(8)=0;
  lddomnode193(8)=130.61;
  lddomnode171(8)=162.98;


  
  tstr193{9}='12_234057';
tsec193(9)=82850;
 tstr171{9}='12_234423';
tsec171(9)=83676;
  lheight193(9)=108.0;
  lheight171(9)=121.3;
  llength193(9)=106.0;
  llength171(9)=101.49;
  lwidth193(9)=5.3;
  lwidth171(9)=6.6;
% lintensity193(9)=0;
% lintensity171(9)=0;
  lddomnode193(9)=134.08;
  lddomnode171(9)=157.0;
  
  


 tstr193{10}='13_000843';
tsec193(10)=84516;
 tstr171{10}='13_001223';
tsec171(10)=85356;
  lheight193(10)=115.48;
  lheight171(10)=117.44;
  llength193(10)=107.0;
  llength171(10)=95.34;
  lwidth193(10)=5.24;
  lwidth171(10)=4.725;
% lintensity193(10)=0;
% lintensity171(10)=0;
  lddomnode193(10)=149.59;
  lddomnode171(10)=150.42;

tstr193{11}='13_002431';
tsec193(11)=85464;
tstr171{11}='13_002723';
tsec171(11)=86256;
  lheight193(11)=94.51;
  lheight171(11)=125.4;
  llength193(11)=132.5;
  llength171(11)=108.5;
  lwidth193(11)=3.7;
  lwidth171(11)=5.3;
% lintensity193(11)=0;
% lintensity171(11)=0;
  lddomnode193(11)=114.28;
  lddomnode171(11)=156.77;
  
  tstr193{12}='';
tsec193(12)=0;
 tstr171{12}='';
tsec171(12)=0;
lheight193(12)=149.09;
lheight171(12)=123.22;
llength193(12)=127.19;
llength171(12)=116.18;
lwidth193(12)=2.45;
lwidth171(12)=3.58;
% % lintensity193(12)=0;
% % lintensity171(12)=0;
% lddomnode193(12)=;
% lddomnode171(12)=;

tstr193{13}='';
tsec193(13)=0;
tstr171{13}='';
tsec171(13)=0;
lheight193(13)=114.33;
lheight171(13)=116.74;
llength193(13)=95.08;
llength171(13)=95.99;
lwidth193(13)=6.02;
lwidth171(13)=5.0;
% lintensity193(13)=0;
% lintensity171(13)=0;
% lddomnode193(13)=;
% lddomnode171(13)=;