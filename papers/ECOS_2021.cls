%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coupled class for LaTeX files                 v1.1 Feb-23-2000 %%
%%                                                                %%
%% The regular article class commands still exist and may be used.%%
%% 10pt, 11pt and twocolumn options are disabled.                 %%
%%                                                                %%
%% \maketitle make the first page automatically                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}

\DeclareOption{10pt}{\OptionNotUsed}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption{twocolumn}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions\relax

\LoadClass[a4paper,12pt]{article}

\RequirePackage{calc}
\RequirePackage{indentfirst}
\RequirePackage{enumitem}  % added, required to adjust vertical list spacing

\setlength{\oddsidemargin}{-5mm}
\setlength{\evensidemargin}{-5mm}
\setlength{\voffset}{-2pt}  % added
\setlength{\topmargin}{0mm}  % added
\setlength{\headheight}{0mm}  % added
\setlength{\headsep}{0mm}  % added
\setlength{\textheight}{257mm}  % adjusted, was 220
\setlength{\textwidth}{170mm}
\setlength{\textfloatsep}{9pt plus 2.0pt minus 4.0pt}  % adjusted, set vertical spacing between text and floats
\setlength{\intextsep}{9pt plus 2.0pt minus 2.0pt}  % added, set vertical distance between text and float
\setlength{\abovedisplayskip}{5pt}
\setlength{\belowdisplayskip}{5pt}
\setlength{\parindent}{0mm}
\setlength{\abovecaptionskip}{3pt}  % adjusted, set vertical spacing above captions
\setlength{\belowcaptionskip}{3pt}  % adjusted, set vertical spacing below captions
\setlength{\itemsep}{0pt}

\setlist{itemsep=1pt, parsep=1pt, topsep=4pt}  % added, set vertical spacing between list items and surrounding paragraphs

%% Title
\renewcommand{\maketitle}{%
 \newpage
% \thispagestyle{myheadings}

 \begin{minipage}[t][1cm][t]{17cm}
  \begin{flushright}
   \scriptsize
\vspace{-34pt}  % adjusted, was -70pt
PROCEEDINGS OF ECOS 2021 - THE $34^ \mathrm{th}$ INTERNATIONAL CONFERENCE ON \\
EFFICIENCY, COST, OPTIMIZATION, SIMULATION AND ENVIRONMENTAL IMPACT OF ENERGY SYSTEMS \\
JUNE 28-JULY 2, 2021, TAORMINA, ITALY

 %  \
 \end{flushright}
 \end{minipage}

 \begin{center}%
  \vspace{-13pt}% adjusted, was 30pt, move title upwards to match Word template title position
  {%3
   \fontsize{18}{18}\selectfont
   \textbf{\@title}
   \par
  }%3
  %\vspace{18pt}%
 \vspace{10pt}%
  {%4
\fontsize{12}{14}\selectfont
   %\bfseries\it\@author
   \textbf{\textit\@author}
   \par
  }%4
  %\vspace{14pt}%
  \vspace{16pt}%
  {%5
\fontsize{10}{12}\selectfont
   %\def\and{\\\vspace{14pt}}
\def\and{\\\vspace{2pt}}
   %\normalfont\@address
   \normalfont\it\@address
   \par
  }%5
 \end{center}

 %\vspace{18pt}%
 \vspace{8pt}%
 \hspace{0cm}\parbox{17cm}{\textbf{\sffamily Abstract:\\}{\textsf\@abstract}}
 \par
  %\vspace{18pt}%
 \vspace{11pt}%
 \hspace{0cm}\parbox{17cm}{\textbf{\sffamily Keywords:\\}{\textsf\@keywords}}
 \par
\vspace{6pt}
 %\vfill
 %\newpage
}% End of \@maketitle

\def\author{\@ifnextchar[{\@getshortauthor}{\@getauthor}}
\def\@getauthor#1{\gdef\@author{#1}\gdef\@shortauthor{#1}}
\def\@getshortauthor[#1]#2{\gdef\@shortauthor{#1}\gdef\@author{#2}}

\newcommand{\address}[1]{\gdef\@address{#1}}

\newcommand{\keywords}[1]{\gdef\@keywords{#1}}

\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}

\newcommand{\heading}[1]{\gdef\@heading{#1}}

%\newcommand{\author}[1]{\gdef\@author{#1}}

%%%%% Headings definition
%\def\ps@myheadings{%

%    \def\@oddfoot{\normalfont\hfil\thepage}
%    \let\@evenfoot\@oddfoot
%    \def\@evenhead{%
%       \parbox[t][11.5pt][t]{\textwidth}{%
%       \parbox[t][12.5pt][t]{\textwidth}{%
%        \centering
%        \normalfont
%        \fontsize{10}{10}\selectfont
%        \ifodd\value{page}\relax
%            \@title\\
%        \else
%            \@author\\
%        \fi%
%        \vspace{-.2cm}
%        \rule{\textwidth}{0.5pt}%
%         }
%    }%
%    \let\@oddhead\@evenhead
%    \let\@mkboth\@gobbletwo
%    \let\sectionmark\@gobble
%    \let\subsectionmark\@gobble
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\sectionLevel{1}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{6pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\bfseries}

\renewcommand\thesection{\arabic{section}.}

\renewcommand \thesubsection{\arabic{section}.\arabic{subsection}.}

\renewcommand \thesubsubsection{\arabic{section}.\arabic{subsection}.\arabic{subsubsection}.}

\renewcommand{\section}{%
  \@startsection{section}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%

\def\sectionLevel{2}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{6pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\large\bfseries}

\def\sectionLevel{3}
\def\sectionIndent{0pt}
\def\sectionSpaceBefore{6pt}
\def\sectionSpaceAfter{6pt}
\def\sectionStyle{\bfseries }

\renewcommand{\subsection}{%
  \@startsection{subsection}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%

\renewcommand{\subsubsection}{%
  \@startsection{subsubsection}{\sectionLevel}{\sectionIndent}{\sectionSpaceBefore}%
                {\sectionSpaceAfter}{\sectionStyle}}%

\renewcommand{\@makecaption}[2]{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\normalsize\fontsize{10}{12}\selectfont {\bf #1}: #2}%
  \ifdim \wd\@tempboxa >\hsize
    \normalsize\fontsize{10}{12}\selectfont {\bf #1}: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}


\renewenvironment{thebibliography}[1]
  {%\newpage
 %  \section*{\MakeUppercase\refname\@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
 \section*{\refname\@mkboth{\refname}{\refname}}%
   \list{\@biblabel{\@arabic\c@enumiv}}%
        {\settowidth\labelwidth{\@biblabel{#1}}%
         \leftmargin\labelwidth
         \advance\leftmargin\labelsep
         \@openbib@code
         \usecounter{enumiv}%
         \let\p@enumiv\@empty
         \renewcommand\theenumiv{\@arabic\c@enumiv}}%
   \setlength{\itemsep}{2pt}  % added, sets vertical space between bib items to match Word template
   \setlength{\parskip}{1pt}  % added, sets vertical space between bib items/paragraph to match Word template
   \sloppy
   \clubpenalty4000
   \@clubpenalty \clubpenalty
   \widowpenalty4000%
   \sfcode`\.\@m}
  {\def\@noitemerr
   {\@latex@warning{Empty `thebibliography' environment}}%
   \endlist}
   \topmargin-.4cm

\renewcommand{\descriptionlabel}[1]{\hspace{1cm}\textit{#1}}

%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\labelitemii}{$\circ$}
%\renewcommand{\labelitemii}{$\square$}
\let\tempone\itemize
\let\temptwo\enditemize
\renewenvironment{itemize}{\tempone\addtolength{\itemsep}{0\baselineskip}}{\temptwo}

%%%%%%%%%%%%%%%%%%%%%
\newenvironment{nomenclature}{\itemize\addtolength{\itemsep}{0.5pt}}{\enditemize}

%\newenvironment{mydescription}{%
%   \renewcommand\descriptionlabel{\h\hspace{1cm}}
%    \renewcommand\descriptionlabel{\hspace{1cm}}
%   \renewcommand\leftmargin{1pt}
%   \begin{description}%
%}{%
%   \end{description}%
%}

\newenvironment{mydescription}[1]
  {\begin{list}{}%
   {\renewcommand\makelabel[1]{##1:\hfill}%
   \settowidth\labelwidth{\makelabel{#1}}%
   \setlength\leftmargin{\labelwidth}
   \addtolength\leftmargin{\labelsep}}}
  {\end{list}}

\AtBeginDocument{%
%  \pagestyle{myheadings}
  \maketitle
  \let\maketitle\relax
%\pagestyle{empty}
  \thispagestyle{plain}
  \pagenumbering{arabic}
}

\endinput
