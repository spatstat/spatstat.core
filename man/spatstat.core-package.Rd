\name{spatstat.core-package}
\alias{spatstat.core-package} 
\alias{spatstat.core} 
\docType{package}
\title{The spatstat.core Package}
\description{
  The \pkg{spatstat.core} package
  belongs to the \pkg{spatstat} family of packages.
  It contains the core functionality
  for statistical analysis and modelling of spatial data.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.


  This sub-package \pkg{spatstat.core} contains all the main
  user-level functions
  that perform statistical analysis and modelling of spatial data.

  (The main exception is that functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.)
}
\section{Structure of the spatstat family}{  
  The orginal \pkg{spatstat} package grew to be very large.
  It has now been divided into several \bold{sub-packages}:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.core} containing the main functionality
    for statistical analysis and modelling of spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.
  For an overview of all the functions available in these sub-packages,
  see the help file for \pkg{spatstat} in the \pkg{spatstat} package,

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Julian Gilbey, 
  Yongtao Guan,
  Ute Hahn,
  Kassel Hingee,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  For comments, corrections, bug alerts and suggestions, we thank
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

