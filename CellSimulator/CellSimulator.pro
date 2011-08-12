#-------------------------------------------------
#
# Project created by QtCreator 2011-07-04T11:34:26
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = CellSimulator
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH = $$PWD/Includes/

DEPENDPATH = $$PWD/Includes/

SOURCES += \
    Sources/NK.cpp \
    Sources/Media.cpp \
    Sources/main.cpp \
    Sources/LT.cpp \
    Sources/experiment1.cpp \
    Sources/Cell_simulator.cpp \
    Sources/APC.cpp \
    Sources/SimParameters.cpp \
    Sources/Results.cpp \
    Sources/ResultsSimulator.cpp \
    Sources/Measurement.cpp \
    Sources/Experiment.cpp \
    Sources/LevenbergMarquardt.cpp \
    Sources/MatrixInverse.cpp

HEADERS += \
    Includes/NK.h \
    Includes/Media.h \
    Includes/LT.h \
    Includes/experiment1.h \
    Includes/Cell_simulator.h \
    Includes/APC.h \
    Includes/SimParameters.h \
    Includes/Results.h \
    Includes/Measurement.h \
    Includes/ResultsSimulator.h \
    Includes/Treatment.h \
    Includes/Experiment.h \
    Includes/LevenbergMarquardt.h \
    Includes/MatrixInverse.h
