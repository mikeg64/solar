%`jaa.cls' 
% LaTeX2e class file for Journal of Astrophysics and Astronomy
% Copyright Indian Academy of Sciences 1996,2003
%
% Author: Dipankar Bhattacharya, Raman Research Institute
%         Bangalore 560080, India.
%
% jaa.cls          Converted to LaTeX2e class
% release 1.0
% 19 April 2003
%
% jaa.sty
% version 2.3      Modified \@startsection,\@sect,\@seccntformat included
% 20 July 1996     for LaTeX2e compatibilty mode.  Fixed a bug in .aux
%                  file output.  Now detects OFSS/NFSS.  Times font can be
%                  loaded through "times" documentstyle option if 
%                  LaTeX2e installation with "times.sty" is available.
%
% release 2.2      Small adjustments in footnote skip, footnote rule, 
% 03 July 1996     float top and bottom glue.
%
% release 2.1      Side margins modified (01 July '96)
%
% release 2.0      Set in 11 point instead of 10 point style.  Corrected
% 22 June 1996     minor formatting discrepancies.  Times fonts now 
%                  loaded through customized jplain.fmt.  "plate"
%                  option abolished.
%
% release 1.1      Included "plate" option and Times Roman fonts.
% 18 June 1996     Will need jaa_press.tex to be included in the
%                  document to be processed for the Times fonts to
%                  take effect.  Without plate option normal cmr
%                  fonts will be used.  "plate" option should be
%                  used only at the JAA editorial office.
%
% release 1.0      First release of basic `jaa' style,  Including
% 15 June 1996     documentstyle options "doublespacing" and 
%                  "landscape".  Page styles, Abstract and Keywords
%                  environments, maketitle, float fractions, sectioning,
%                  paragraphing, bibliography, page sizes, captions 
%                  and type styles customized.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% article.sty 16-Mar-88
% rri.sty 02-Feb-95 D. Bhattacharya
% last modified 07-Feb-95 DB
%
% Modified to create JAA journal style 
% with material adapted from mn.sty (authors: Sue Brooks, Alison Woollatt 
%                    and Mark Reed, copyright Cambridge University Press)
% D. Bhattacharya 13 June 1996
% First release 15 June 1996
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{jaa}
            [2003/04/19 v1.0 Journal of Astrophysics and Astronomy]

\newif\ifdoublesp@cing
\newif\ifja@landscape
%
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{draft}{\setlength\overfullrule{5pt} \doublesp@cingtrue}
\DeclareOption{final}{\setlength\overfullrule{0pt} \doublesp@cingfalse}
\DeclareOption{doublespacing}{\doublesp@cingtrue}
\DeclareOption{landscape}{\ja@landscapetrue}
\DeclareOption{portrait}{\ja@landscapefalse}
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
%
\ExecuteOptions{a4paper,portrait,twoside,onecolumn,final}
\ProcessOptions
%
%
% Trailing dot added to section numbers, but not subsection etc
%

\newif\if@levelone

\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindenttrue
   \ifnum #2=\@ne  % if level = 1 (section)
     \global\@levelonetrue
   \else
     \global\@levelonefalse
   \fi
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
%
\def\@seccntformat#1{\csname the#1\endcsname\if@levelone.\fi\hskip 1em}
%
\def\@sect#1#2#3#4#5#6[#7]#8{%
    \ifnum #2>\c@secnumdepth
      \def\@svsec{}%
    \else
      \refstepcounter{#1}%
      \let\@@protect\protect
      \def\protect{\noexpand\protect\noexpand}%
      \edef\@svsec{\@seccntformat{#1}}%
      \let\protect\@@protect
    \fi
    \@tempskipa #5\relax
    \ifdim \@tempskipa>\z@
      \begingroup #6\relax
        \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
      \endgroup
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else %
      \protect\numberline{\csname the#1\endcsname}\fi #7}%
    \else
      \def\@svsechd{#6\hskip #3\relax %
        \@svsec #8 \csname #1mark\endcsname {#7}%
        \addcontentsline{toc}{#1}{\ifnum #2>\c@secnumdepth \else %
        \protect\numberline{\csname the#1\endcsname}\fi#7}%
      }\fi%
    \@xsect{#5}}
%
% 
% art11.sty
% rri11.sty 02-Feb-1995 D. Bhattacharya
% Modified and adopted for JAA journal style
%
\typeout{11 point type style}

\lineskip 1pt \normallineskip 1pt
%
\ifdoublesp@cing
 \def\baselinestretch{1.7}
\else
  \def\baselinestretch{1}
\fi
%
\renewcommand\normalsize{%
   \@setfontsize\normalsize\@xipt{13.2}%
   \abovedisplayskip 11\p@ \@plus3\p@ \@minus6\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6.5\p@ \@plus3.5\p@ \@minus3\p@
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\normalsize
\newcommand\small{%
   \@setfontsize\small\@xpt\@xiipt
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 6\p@ \@plus2\p@ \@minus2\p@
               \parsep 3\p@ \@plus2\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\footnotesize{%
   \@setfontsize\footnotesize\@ixpt{11}%
   \abovedisplayskip 8\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\newcommand\scriptsize{\@setfontsize\scriptsize\@viiipt{9.5}}
\newcommand\tiny{\@setfontsize\tiny\@vipt\@viipt}
\newcommand\large{\@setfontsize\large\@xiipt{14}}
\newcommand\Large{\@setfontsize\Large\@xivpt{18}}
\newcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\newcommand\huge{\@setfontsize\huge\@xxpt{25}}
\newcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
%
% Margins modified
%
\oddsidemargin 47pt \evensidemargin 47pt
\marginparwidth 90pt 
\marginparsep 10pt 
%
% Sizes modified
%
\ifja@landscape
  \textheight 361pt \textwidth 567.6pt
\else
  \textheight 567.6pt \textwidth 361pt
\fi
%
 \topmargin 27pt \headheight 12pt \headsep 13pt \footskip 17.5pt
\advance\textheight by \topskip
\columnsep 10pt \columnseprule 0pt 
%
% Footnoteskip modified
%
\footnotesep 7.7pt 
\skip\footins 20pt plus 12pt minus 1pt 
%
\floatsep 12pt plus 2pt minus 2pt \textfloatsep 20pt plus 2pt minus 4pt
\intextsep 12pt plus 2pt minus 2pt \dblfloatsep 12pt plus 2pt
minus 2pt \dbltextfloatsep 20pt plus 2pt minus 4pt  
\setlength\maxdepth{.5\topskip}
%
% float fills modified from 1fil to 0fil for top and bottom
%
\@fptop 0pt plus 0fil \@fpsep 8pt plus 2fil \@fpbot 0pt plus 0fil 
\@dblfptop 0pt plus 0fil \@dblfpsep 8pt plus 2fil \@dblfpbot 0pt plus 0fil
\marginparpush 5pt 
%
% Paragraphing glue and penalty modified
%
\parskip \z@ plus .1pt \parindent 17pt \partopsep 0pt plus .1pt 
\@lowpenalty 51 \@medpenalty 151 \@highpenalty 301 
\@beginparpenalty -\@lowpenalty \@endparpenalty -\@lowpenalty 
\@itempenalty -\@lowpenalty 
\clubpenalty=0 \widowpenalty=10000
%

\def\part{\par \addvspace{4ex} \@afterindentfalse \secdef\@part\@spart} 
\def\@part[#1]#2{\ifnum \c@secnumdepth >\m@ne \refstepcounter{part}
\addcontentsline{toc}{part}{\thepart \hspace{1em}#1}\else
\addcontentsline{toc}{part}{#1}\fi { \parindent 0pt \raggedright 
 \ifnum \c@secnumdepth >\m@ne \Large \bf Part \thepart \par\nobreak \fi \huge
\bf #2\markboth{}{}\par } \nobreak \vskip 3ex \@afterheading } 
\def\@spart#1{{\parindent 0pt \raggedright 
 \huge \bf 
 #1\par} \nobreak \vskip 3ex \@afterheading } 
%

\setcounter{secnumdepth}{3}

\def\appendix{\par
 \setcounter{section}{0}
 \setcounter{subsection}{0}
 \def\thesection{\Alph{section}}}

\leftmargini 2.5em
\leftmarginii 2.2em \leftmarginiii 1.87em \leftmarginiv 1.7em \leftmarginv 1em
\leftmarginvi 1em
\leftmargin\leftmargini
\labelsep .5em
\labelwidth\leftmargini\advance\labelwidth-\labelsep
\def\@listI{\leftmargin\leftmargini \parsep 4.5pt plus 2pt minus 1pt\topsep
9pt plus 3pt minus 5pt\itemsep 4.5pt plus 2pt minus 1pt}
\let\@listi\@listI
\@listi 
\def\@listii{\leftmargin\leftmarginii
 \labelwidth\leftmarginii\advance\labelwidth-\labelsep
 \topsep 4.5pt plus 2pt minus 1pt
 \parsep 2pt plus 1pt minus 1pt
 \itemsep \parsep}
\def\@listiii{\leftmargin\leftmarginiii
 \labelwidth\leftmarginiii\advance\labelwidth-\labelsep
 \topsep 2pt plus 1pt minus 1pt 
 \parsep \z@ \partopsep 1pt plus 0pt minus 1pt
 \itemsep \topsep}
\def\@listiv{\leftmargin\leftmarginiv
 \labelwidth\leftmarginiv\advance\labelwidth-\labelsep}
\def\@listv{\leftmargin\leftmarginv
 \labelwidth\leftmarginv\advance\labelwidth-\labelsep}
\def\@listvi{\leftmargin\leftmarginvi
 \labelwidth\leftmarginvi\advance\labelwidth-\labelsep}
%
% end of rri11.sty
%
%
\relax
\def\labelenumi{\theenumi.} 
\def\theenumi{\arabic{enumi}} 
\def\labelenumii{(\theenumii)}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv.}
\def\theenumiv{\Alph{enumiv}} 
\def\p@enumiv{\p@enumiii\theenumiii}

\def\labelitemi{$\bullet$}
\def\labelitemii{\bf --}
\def\labelitemiii{$\ast$}
\def\labelitemiv{$\cdot$}

\def\verse{\let\\=\@centercr 
 \list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent 
 \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse\endlist
\def\quotation{\list{}{\listparindent 1.5em
 \itemindent\listparindent
 \rightmargin\leftmargin \parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist
\def\quote{\list{}{\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

\def\descriptionlabel#1{\hspace\labelsep \bf #1}
\def\description{\list{}{\labelwidth\z@ \itemindent-\leftmargin
 \let\makelabel\descriptionlabel}}
\let\enddescription\endlist
%
\def\theequation{\arabic{equation}}

\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
 \else \newpage \fi \thispagestyle{empty}\c@page\z@}
\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

\arraycolsep 5pt \tabcolsep 6pt \arrayrulewidth .4pt \doublerulesep 2pt 
\tabbingsep \labelsep 

\skip\@mpfootins = \skip\footins
\fboxsep = 3pt \fboxrule = .4pt 
%
\newcounter{part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

\def\thepart{\Roman{part}} 
%
\def\thesection {\arabic{section}} 
\def\thesubsection {\mbox{\rm \thesection.\arabic{subsection}}} 
\def\thesubsubsection {\mbox{\rm \thesubsection .\arabic{subsubsection}}}
%
\def\theparagraph {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph {\theparagraph.\arabic{subparagraph}}
%
% Sectioning style modified
%
\def\section{\@startsection{section}{1}{\z@}%
 {-3.5ex plus -1ex minus -.2ex}{2.3ex plus .2ex}{\reset@font\large\bf}}
\def\subsection{\@startsection{subsection}{2}{\z@}%
 {-3.25ex plus -1ex minus -.2ex}{1.5ex plus .2ex}{\reset@font\normalsize\it}}
\def\subsubsection{\@startsection{subsubsection}{3}{\z@}%
 {-3.25ex plus -1ex minus -.2ex}{1.5ex plus .2ex}{\reset@font\normalsize\it}}
\def\paragraph{\@startsection{paragraph}{4}{\z@}{3.25ex plus 1ex minus .2ex}%
 {-1em}{\reset@font\normalsize\bf}}
\def\subparagraph{\@startsection{subparagraph}{4}{\parindent}%
 {3.25ex plus 1ex minus .2ex}{-1em}{\reset@font\normalsize\it}}
%
\def\@pnumwidth{1.55em}
\def\@tocrmarg {2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{3}

\def\tableofcontents{\section*{Contents\@mkboth{CONTENTS}{CONTENTS}}
 \@starttoc{toc}}
\def\l@part#1#2{\addpenalty{\@secpenalty}
 \addvspace{2.25em plus 1pt} \begingroup
 \@tempdima 3em \parindent \z@ \rightskip \@pnumwidth \parfillskip
-\@pnumwidth 
 {\large \bf \leavevmode #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
 \nobreak \endgroup}
\def\l@section#1#2{\addpenalty{\@secpenalty} \addvspace{1.0em plus 1pt}
\@tempdima 1.5em \begingroup
 \parindent \z@ \rightskip \@pnumwidth 
 \parfillskip -\@pnumwidth 
 \bf \leavevmode \advance\leftskip\@tempdima \hskip -\leftskip #1\nobreak\hfil
\nobreak\hbox to\@pnumwidth{\hss #2}\par
 \endgroup}
 
\def\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\def\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\def\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\def\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\def\listoffigures{\section*{List of Figures\@mkboth
 {LIST OF FIGURES}{LIST OF FIGURES}}\@starttoc{lof}}
\def\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\def\listoftables{\section*{List of Tables\@mkboth
 {LIST OF TABLES}{LIST OF TABLES}}\@starttoc{lot}}
\let\l@table\l@figure
%
%  thebibliography definition modified
%
%
\newcounter{dummy}
\newlength{\leftadvance}\setlength{\leftadvance}{3em}
\def\thebibliography#1{\section*{References}\list
 {\mbox{}}{\usecounter{dummy}
 \leftmargin\leftadvance \itemindent-\leftadvance
 \labelwidth\z@ \itemsep\z@}
 \reset@font\small
%\small 
 \parindent\z@ \parskip\z@ plus .1pt \relax
 \def\newblock{\hskip .11em plus .33em minus .07em}
 \sloppy\clubpenalty4000\widowpenalty4000
 \sfcode`\.=1000\relax}
\let\endthebibliography=\endlist
%
\newif\if@restonecol
\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
\columnseprule \z@
\columnsep 35pt\twocolumn[\section*{Index}]
 \@mkboth{INDEX}{INDEX}\thispagestyle{plain}\parindent\z@
 \parskip\z@ plus .3pt\relax\let\item\@idxitem}
\def\@idxitem{\par\hangindent 40pt}
\def\subitem{\par\hangindent 40pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}
%
\def\footnoterule{\kern-3\p@ 
 \hrule width 57pt
 \kern 2.6\p@} 
\long\def\@makefntext#1{\parindent 1em\noindent 
 \hbox to 1.8em{\hss$^{\@thefnmark}$}#1}
%
% Float fractions modified
%
\setcounter{topnumber}{2}
\def\topfraction{.9}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.9}
\setcounter{totalnumber}{3}
\def\textfraction{.1}
\def\floatpagefraction{.9}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.9}
\def\dblfloatpagefraction{.9}
%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
%
% font size in figure caption changed to small -- DB 7/2/95
%
\def\@captionfont{\small}
%
% Makecaption customized
%
\long\def\@makecaption#1#2{
 \vskip 10pt 
 \setbox\@tempboxa\hbox{\@captionfont {\bf #1.} #2}
 \ifdim \wd\@tempboxa >\hsize {\@captionfont {\bf #1.} #2}\par \else \hbox
to\hsize{\hfil\box\@tempboxa\hfil} 
 \fi}
%
\newcounter{figure}
\def\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{{Figure \thefigure}}
\def\figure{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}

\newcounter{table}
\def\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{{Table \thetable}}
\def\table{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}
%
\def\thanks#1{\footnotemark\hbox{~}\begingroup
\def\protect{\noexpand\protect\noexpand}\xdef\@thanks{\@thanks
 \protect\footnotetext[\the\c@footnote]{#1}}\endgroup}
\def\@thanks{}
%
\def\maketitle{\par
 \begingroup
 \def\thefootnote{\fnsymbol{footnote}}
 \def\@makefnmark{\hbox
 to 0pt{$^{\@thefnmark}$\hss}}
 \if@twocolumn
 \twocolumn[\@maketitle]
 \else
  \newpage
  \global\@topnum\z@
  \@maketitle
 \fi
 \thispagestyle{jaafirst}
 \@thanks
 \endgroup
 \setcounter{footnote}{0}
 \let\maketitle\relax
 \let\@maketitle\relax
 \gdef\@author{}
 \gdef\@title{}
 \gdef\@thanks{}
 \let\thanks\relax
}

\def\newauthor{%
  \end{author@tabular}\par
  \begin{author@tabular}[t]{@{}l@{}}%
}
\def\and{%
  \end{author@tabular}\vskip 3pt\par
 \begin{author@tabular}[t]{@{}l@{}}%
}
\def\@maketitle{\newpage
 \vspace*{20pt}
 {\raggedright \sloppy
  {\reset@font\Large\bf \@title \par}
% {\Large\bf \@title \par}
  \vskip 7pt
  {\reset@font\Large
%  {\Large
   \begin{author@tabular}[t]{@{}l@{}}\@author
   \end{author@tabular}\par}
  \vskip 7pt
 }
 \par\noindent
 {\reset@font\small \@date \par}
% {\small \@date \par}
 \vskip 6pt
}
%
\def\abstract{\if@twocolumn
\section*{Abstract}
\else
\quotation {\noindent \bf Abstract.\ }
 \rm
\fi}
\def\endabstract{\if@twocolumn\else\endquotation\fi}
%
\def\keywords{\if@twocolumn
\section*{Keywords}
\else
\quotation {\noindent {\it Key words}:}
 \rm
\fi}
\def\endkeywords{\if@twocolumn\else\endquotation\vspace{3mm}\fi}
%
\def\author@tabular{\def\@halignto{}\@authortable}
\let\endauthor@tabular=\endtabular
%
\def\author@tabcrone{\vspace{2pt}{\ifnum0=`}\fi%
  \@xtabularcr[-7pt]\reset@font\small\it
% \@xtabularcr[-7pt]\small\it
  \let\\=\author@tabcrtwo\ignorespaces}
%
\def\author@tabcrtwo{{\ifnum0=`}\fi\@xtabularcr[-7pt]\reset@font\small\it
%\def\author@tabcrtwo{{\ifnum0=`}\fi\@xtabularcr[-7pt]\small\it
 \let\\=\author@tabcrtwo\ignorespaces}
%
\def\@authortable{\leavevmode \hbox \bgroup $\let\@acol\@tabacol
 \let\@classz\@tabclassz \let\@classiv\@tabclassiv
 \let\\=\author@tabcrone\ignorespaces \@tabarray}
%
\mark{{}{}} 

\if@twoside \def\ps@headings{\let\@mkboth\markboth
\def\@oddfoot{}\def\@evenfoot{}\def\@evenhead{\rm \thepage\hfil \sl
\leftmark}\def\@oddhead{\hbox{}\sl \rightmark \hfil
\rm\thepage}\def\sectionmark##1{\markboth {\uppercase{\ifnum \c@secnumdepth
>\z@
 \thesection\hskip 1em\relax \fi ##1}}{}}\def\subsectionmark##1{\markright
{\ifnum \c@secnumdepth >\@ne
 \thesubsection\hskip 1em\relax \fi ##1}}}
\else \def\ps@headings{\let\@mkboth\markboth
\def\@oddfoot{}\def\@evenfoot{}\def\@oddhead{\hbox {}\sl \rightmark \hfil
\rm\thepage}\def\sectionmark##1{\markright {\uppercase{\ifnum \c@secnumdepth
>\z@
 \thesection\hskip 1em\relax \fi ##1}}}}
\fi
\def\ps@myheadings{\let\@mkboth\@gobbletwo
\def\@oddhead{\hbox{}\sl\rightmark \hfil
\rm\thepage}\def\@oddfoot{}\def\@evenhead{\rm \thepage\hfil\sl\leftmark\hbox
{}}\def\@evenfoot{}\def\sectionmark##1{}\def\subsectionmark##1{}}
%
% JAA page style for pages other than first page of article
%
\def\ps@jaaheadings{\let\@mkboth\@gobbletwo
 \def\@oddhead{\hbox{}\hfill{\it\@shorttitle}\hfill\hbox{}\rm\thepage}
 \def\@oddfoot{}
 \def\@evenhead{\rm\thepage\hbox{}\hfill{\it\@shortauthor}\hfill\hbox{}}
 \def\@evenfoot{}
 \def\sectionmark##1{\markboth{##1}{}}
 \def\subsectionmark##1{\markright{##1}}}
%
% JAA first page style
%
\gdef\@author{\mbox{}}
\gdef\@shortauthor{}
\def\author{\@ifnextchar [{\@authortwo}{\@authorone}}
\def\@authortwo[#1]#2{\gdef\@author{#2}\gdef\@shortauthor{#1}}
\def\@authorone#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\gdef\@title{\mbox{}}
\gdef\@shorttitle{}
\def\title{\@ifnextchar [{\@titletwo}{\@titleone}}
\def\@titletwo[#1]#2{\gdef\@title{#2}\gdef\@shorttitle{#1}}
\def\@titleone#1{\gdef\@title{#1}\gdef\@shorttitle{#1}}
%
\def\volume#1{\gdef\@volume{#1}}
\gdef\@volume{00}
\def\pubyear#1{\gdef\@pubyear{#1}}
\gdef\@pubyear{0000}
\def\pagerange#1{\gdef\@pagerange{#1}}
\gdef\@pagerange{000--000}
\def\journal#1{\gdef\@journal{#1}}
\gdef\@journal{{J.\ Astrophys.\ Astr.} (\@pubyear) {\bf \@volume}, \@pagerange}
%
\def\ps@jaafirst{\let\@mkboth\@gobbletwo
    \def\@oddhead{\small\@journal\hfil}
    \def\@oddfoot{\hfil \rm\thepage}
    \def\@evenhead{\small\@journal\hfil}
    \def\@evenfoot{\hfil \rm\thepage}}
%
\def\today{\number\year\ \ifcase\month\or
 January\or February\or March\or April\or May\or June\or
 July\or August\or September\or October\or November\or December
 \fi \ \number\day}
%
\ps@jaaheadings
\pagenumbering{arabic} \onecolumn \if@twoside\else\raggedbottom\fi 
\endinput
%%
%% End of file `jaa.cls'
