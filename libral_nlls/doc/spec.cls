%%
%% This is file `hslspec.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% classes.dtx  (with options: `hslspec')
%% 
%% This is a generated file.
%% 
%% Copyright 1993 1994 1995 1996 1997 1998 1999
%% The LaTeX3 Project and any individual authors listed elsewhere
%% in this file.
%% 
%% This file is part of the LaTeX2e system.
%% ----------------------------------------
%% 
%% It may be distributed under the terms of the LaTeX Project Public
%% License, as described in lppl.txt in the base LaTeX distribution.
%% Either version 1.0 or, at your option, any later version.
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{hslspec}
              [1999/01/07 v1.4a
 Standard LaTeX document class]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\fi
\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
%\ExecuteOptions{letterpaper,10pt,oneside,onecolumn,final}
\ExecuteOptions{a4paper,10pt,twoside,onecolumn,final}
\ProcessOptions
\input{size1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{2\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.5}
\if@twoside
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\slshape\leftmark}%
      \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\MakeUppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\quad
        \fi
        ##1}}}
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markright {\MakeUppercase{%
        \ifnum \c@secnumdepth >\m@ne
          \thesection\quad
        \fi
        ##1}}}}
\fi
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }
  \if@titlepage
  \newcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \null\vfil
  \vskip 60\p@
  \begin{center}%
    {\LARGE \@title \par}%
    \vskip 3em%
    {\large
     \lineskip .75em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
      \vskip 1.5em%
    {\large \@date \par}%       % Set date in \large size.
  \end{center}\par
  \@thanks
  \vfil\null
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\else
\newcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\LARGE \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    \vskip 1em%
    {\large \@date}%
  \end{center}%
  \par
  \vskip 1.5em}
\fi
\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand \thepart {\@Roman\c@part}
\renewcommand \thesection {\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\newcommand\part{%
   \if@noskipsec \leavevmode \fi
   \par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname~\thepart
       \par\nobreak
     \fi
     \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\normalfont\large\bfseries}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \listparindent1.0em
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\if@titlepage
  \newenvironment{abstract}{%
      \titlepage
      \null\vfil
      \@beginparpenalty\@lowpenalty
      \begin{center}%
        \bfseries \abstractname
        \@endparpenalty\@M
      \end{center}}%
     {\par\vfil\null\endtitlepage}
\else
  \newenvironment{abstract}{%
      \if@twocolumn
        \section*{\abstractname}%
      \else
        \small
        \begin{center}%
          {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}%
        \end{center}%
        \quotation
      \fi}
      {\if@twocolumn\else\endquotation\fi}
\fi
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\if@compatibility
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\z@
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
    }
\else
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
     \if@twoside\else
        \setcounter{page}\@ne
     \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \gdef\thesection{\@Alph\c@section}}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand \theequation {\@arabic\c@equation}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename
      \@mkboth{%
          \MakeUppercase\listtablename}{\MakeUppercase\listtablename}}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
     {\section*{\refname
        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{References}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi

\newcommand{\hsl}{}

\usepackage{graphicx,fancyhdr,pslatex}

% PDF details
\usepackage{ifpdf}
\usepackage{url}


% Page sizes

\setlength\textheight {21.1cm}
%\setlength\textheight {17.5cm}
\setlength\textwidth {16.5 cm}
\setlength\oddsidemargin {-1.0mm}
\setlength\evensidemargin{-2.0mm}
\setlength\topmargin {14.0mm}
\setlength\headheight {1.05cm}
\setlength\headsep {3.0mm}
\setlength\topskip {5.0mm}
%\setlength\footskip {1.1cm}


% Page style

\pagestyle{fancyplain}
\renewcommand{\headrulewidth}{0.5mm}
\renewcommand{\footrulewidth}{0.5mm}

\newcommand{\licence}{\textbf{\large Licence: 2-clause BSD.}}
\newcommand{\website}{\small \url{http://ccpforge.cse.rl.ac.uk/gf/project/ral_nlls/}}

\fancyhead{} % clear all header fields
\fancyhead[LE,RO]{\Large \textbf{\fullpackagename}}
\fancyhead[LO,RE]{\Large \textbf{\releasename}}

\fancyfoot{} % clear all footer fields
\fancyfoot[RE,RO]{\texttt{\fullpackagename} v\version \\ Documentation date: \today}
\fancyfoot[LE,LO]{\licence\\\website}
\fancyfoot[CE,CO]{\ \\\thepage}

% Redefine the plain style so there is no header, only a footer
\fancypagestyle{plain}{
   \fancyhead{} % clear header
   \renewcommand{\headrulewidth}{0.0mm}
}

%     New commands

\renewcommand{\thesection}{\arabic{section}}
\renewcommand{\theequation}{\arabic{section}.\arabic{equation}}
\newcommand{\releasename}{}

%  Specifics

\newcommand{\itt}[1]{\item[\tt #1]}
\newcommand{\ittf}[1]{\item[\tt \parbox{7mm}{{#1}}]}
\setcounter{secnumdepth}{3}

\newcommand{\subroutine}[1]{\smallskip{\setlength{\labelsep}{0cm}
\setlength{\leftmargini}{2cm}
\begin{itemize}\item[{\tt CAL}]{\tt L #1}\end{itemize}}\smallskip}

\newcommand{\function}[1]{\smallskip{\setlength{\labelsep}{0cm}
\setlength{\leftmargini}{1.5cm}
\begin{itemize}\item[]{\tt #1}\end{itemize}}\smallskip}

% front page header 

\vspace*{-30mm}
\newcommand{\hslheader}{
\ifpdf
\pdfinfo{
   /Author (Science and Technology Facilities Council)
   /Title  (\libraryname_\packagename)
}
\fi
\thispagestyle{plain}

\noindent
\begin{figure}[th]
\includegraphics[height=30pt,trim=0.0cm 0.35cm 0.0cm 0.1cm]{stfc}
\hfill
{\huge \bf \fullpackagename}
\end{figure}
\vspace*{-5.0mm}

\noindent
\rule{\textwidth}{0.5mm}

\noindent
{PACKAGE SPECIFICATION \hfill \releasename}
\vspace*{-1.0mm}

\noindent
\rule{\textwidth}{0.5mm}
\vspace*{-5.0mm}
}

% section names

\newcommand{\resetcounters}{\setcounter{equation}{0} \setcounter{figure}{0} 
 \setcounter{table}{0}}
\newcommand{\hslsummary}{\section{SUMMARY\label{hslsummary}} \resetcounters}
\newcommand{\hslhowto}{\section{HOW TO USE THE PACKAGE\label{hslhowto}}
                        \resetcounters}
\newcommand{\hslgeneral}{\section{GENERAL INFORMATION\label{hslgeneral}}
                         \resetcounters \begin{description}}
\newcommand{\hslmethod}{\end{description}\section{METHOD\label{hslmethod}}
                        \resetcounters}
\newcommand{\hslexample}{\section{EXAMPLE OF USE\label{hslexample}}
                        \resetcounters}

% subsection names

\newcommand{\hslmatrix}{\subsection{Matrix storage formats\label{hslmatrix}}}
\newcommand{\hsltype}{\subsection{The derived data type\label{hsltype}}}
\newcommand{\hsltypes}{\subsection{The derived data types\label{hsltypes}}}
\newcommand{\hslarguments}{\subsection{Argument lists and calling sequences\label{hslarguments}}}
\newcommand{\hslreverse}{\subsection{Reverse communication\label{hslreverse}}}
\newcommand{\hslerrors}{\subsection{Warning and error messages\label{hslerrors}}}
\newcommand{\hslinfo}{\subsection{Information printed\label{hslinfo}}}

% attributes

\newcommand{\hslattributes}{\medskip\noindent{\bf ATTRIBUTES --- }}
\newcommand{\hslversions}[1]{{\bf Version:} #1}
\newcommand{\hslIRDCZ}{{\bf Types: }}
\newcommand{\hsluses}{{\bf Uses: }}
\newcommand{\hslcalls}{{\bf Calls: }}
\newcommand{\hsldate}{{\bf Date: }}
\newcommand{\hsloriginaldate}{{\bf Original date: }}
\newcommand{\hslorigin}{{\bf Origin: }}
\newcommand{\hsllanguage}{{\bf Language: }}
\newcommand{\hslparallelism}{{\bf Parallelism: }}
\newcommand{\hslremark}{{\bf Remark: }}
\newcommand{\hslinterfaces}{{\bf Interfaces: }}

% description names

\newcommand{\hslcommon}{\item[Use of common:]}
\newcommand{\hslworkspace}{\item[Workspace:]}
\newcommand{\hslroutines}{\item[Other routines called directly:]}
\newcommand{\hslmodules}{\item[Other modules used directly:]}
\newcommand{\hslio}{\item[Input/output:]}
\newcommand{\hslrestrictions}{\item[Restrictions:]}
\newcommand{\hslportability}{\item[Portability:]}
\newcommand{\hslchanges}{\item[Changes from Version]}
\newcommand{\hslreference}{\vspace*{5mm} \noindent {\large {\bf Reference:} }}
\newcommand{\hslreferences}{\vspace*{5mm} \noindent {\large {\bf References:} }}

% data types

\newcommand{\integer}{{\tt INTEGER}}
\newcommand{\logical}{{\tt LOGICAL}}
\newcommand{\character}{{\tt CHARACTER}}
\newcommand{\real}{{\tt REAL}}
\newcommand{\realdp}{{\tt REAL} (double precision 
 in {\tt \fullpackagename\_double})}
\newcommand{\complex}{{\tt COMPLEX}}
\newcommand{\complexdp}{{\tt COMPLEX} (double precision complex 
 in {\tt \fullpackagename\_double})}
\newcommand{\optional}{{\tt OPTIONAL}}
\newcommand{\intentin}{{\tt INTENT(IN)}}
\newcommand{\intentinout}{{\tt INTENT(INOUT)}}
\newcommand{\intentout}{{\tt INTENT(OUT)}}
\newcommand{\true}{{\tt .TRUE.}}
\newcommand{\false}{{\tt .FALSE.}}
\newcommand{\restrictions}{{\bf Restrictions: }}
\newcommand{\restriction}{{\bf Restriction: }}

% Other useful definitions

\newcommand{\eqn}[2]{\begin{equation}\label{#1}{#2}\end{equation}}
\newcommand{\disp}[1]{\[{#1}\]}

\newcommand{\vect}[1]{\left(\begin{array}{c}#1\end{array}\right)}
\newcommand{\mat}[2]{\left(\begin{array}{#1}#2\end{array}\right)}
\newcommand{\matinv}[2]{\left(\begin{array}{#1}#2\end{array}\right)^{-1}}
\newcommand{\arr}[2]{\begin{array}{#1}#2\end{array}}

\newcommand{\req}[1]{(\ref{#1})}
\newcommand{\rim}[1]{\;\;\; \mbox{#1} \;\;\;}
\newcommand{\tim}[1]{\;\; \mbox{#1} \;\;}

\newcommand{\minin}[1]{ {\renewcommand{\arraystretch}{0.8} 
                     \begin{array}[t]{c} 
                     \mbox{minimize} \vspace*{-1mm} \\ 
                     \mbox{ $\scriptstyle #1 $ } 
                     \end{array} \;} }
\newcommand{\smallRe}
     {\mbox{\raisebox{-0.8pt}{\footnotesize I\hskip -1.5pt R\hskip -0.5pt}}}
\newcommand{\tinyRe}
     {\mbox{\raisebox{-0.8pt}{\scriptsize I\hskip -1.5pt R\hskip -0.5pt}}}

% Fractions

\newcommand{\sfrac}[2]{{\scriptstyle \frac{#1}{#2}}}
\newcommand{\half}{\sfrac{1}{2}}
\newcommand{\threehalves}{\sfrac{3}{2}}

% Bold mathematics

% syntax: \bm[A-Z,a-z]

\newcommand{\bmA}{{\mathbf{A}}} \newcommand{\bmB}{{\mathbf{B}}}
\newcommand{\bmC}{{\mathbf{C}}} \newcommand{\bmD}{{\mathbf{D}}}
\newcommand{\bmE}{{\mathbf{E}}} \newcommand{\bmF}{{\mathbf{F}}}
\newcommand{\bmG}{{\mathbf{G}}} \newcommand{\bmH}{{\mathbf{H}}}
\newcommand{\bmI}{{\mathbf{I}}} \newcommand{\bmJ}{{\mathbf{J}}}
\newcommand{\bmK}{{\mathbf{K}}} \newcommand{\bmL}{{\mathbf{L}}}
\newcommand{\bmM}{{\mathbf{M}}} \newcommand{\bmN}{{\mathbf{N}}}
\newcommand{\bmO}{{\mathbf{O}}} \newcommand{\bmP}{{\mathbf{P}}}
\newcommand{\bmQ}{{\mathbf{Q}}} \newcommand{\bmR}{{\mathbf{R}}}
\newcommand{\bmS}{{\mathbf{S}}} \newcommand{\bmT}{{\mathbf{T}}}
\newcommand{\bmU}{{\mathbf{U}}} \newcommand{\bmV}{{\mathbf{V}}}
\newcommand{\bmW}{{\mathbf{W}}} \newcommand{\bmX}{{\mathbf{X}}}
\newcommand{\bmY}{{\mathbf{Y}}} \newcommand{\bmZ}{{\mathbf{Z}}}

\newcommand{\bma}{{\mathbf{a}}} \newcommand{\bmb}{{\mathbf{b}}}
\newcommand{\bmc}{{\mathbf{c}}} \newcommand{\bmd}{{\mathbf{d}}}
\newcommand{\bme}{{\mathbf{e}}} \newcommand{\bmf}{{\mathbf{f}}}
\newcommand{\bmg}{{\mathbf{g}}} \newcommand{\bmh}{{\mathbf{h}}}
\newcommand{\bmi}{{\mathbf{i}}} \newcommand{\bmj}{{\mathbf{j}}}
\newcommand{\bmk}{{\mathbf{k}}} \newcommand{\bml}{{\mathbf{l}}}
\newcommand{\bmm}{{\mathbf{m}}} \newcommand{\bmn}{{\mathbf{n}}}
\newcommand{\bmo}{{\mathbf{o}}} \newcommand{\bmp}{{\mathbf{p}}}
\newcommand{\bmq}{{\mathbf{q}}} \newcommand{\bmr}{{\mathbf{r}}}
\newcommand{\bms}{{\mathbf{s}}} \newcommand{\bmt}{{\mathbf{t}}}
\newcommand{\bmu}{{\mathbf{u}}} \newcommand{\bmv}{{\mathbf{v}}}
\newcommand{\bmw}{{\mathbf{w}}} \newcommand{\bmx}{{\mathbf{x}}}
\newcommand{\bmy}{{\mathbf{y}}} \newcommand{\bmz}{{\mathbf{z}}}

\endinput

%%
%% End of file `hslspec.cls'.

