package com.bitbucket.arrc.arrc_acousticroomrecreation;

import android.app.Activity;
import android.content.Intent;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaPlayer;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.os.Handler;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import com.androidplot.xy.LineAndPointFormatter;
import com.androidplot.xy.PanZoom;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.XYSeries;

import java.io.IOException;
import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Arrays;

import org.jtransforms.fft.DoubleFFT_1D;

import static java.lang.Double.NEGATIVE_INFINITY;

public class SweepActivity extends Activity {

    Button mBut;
    Button mRetake;
    XYPlot mPlot;
    TextView mText;

    private static final String LOG_TAG = "ARRC";
    private static final double SPEED_OF_SOUND = 343.0f; // m/s
    private static final int FS = 48000;
    private static final int RECORDER_CHANNELS = AudioFormat.CHANNEL_IN_MONO;
    private static final int RECORDER_AUDIO_ENCODING = AudioFormat.ENCODING_PCM_FLOAT;

    private static final int IR_HSCALE = 100;
    private static final int IR_LEN = FS/20; // 0.05s
    private static final double PEAK_THRESHOLD = 1.06;
    private static final int WINDOW_SIZE = 24; // 200 was good for Dokmanic's 96khz test data
    private static final int WINDOW_POST = FS/40; // 0.025s

    private AudioRecord mRecorder = null;
    private MediaPlayer mPlayer = null;
    private Thread mRecordingThread = null;
    private boolean mIsRecording = false;
    private FloatBuffer mAudioBuf;
    private String mIP;

    private static int NUM_MICS = 5;
    private double[][] mIRs;
    private int[][] mIRpeaks;
    private int[] mIRDirects;
    private int mIRidx = 0;

    private double mDistSpkr;
    private double mHeightDiff;
    private double mXYDiff;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_sweep);
        mBut = findViewById(R.id.sweep_button);
        mRetake = findViewById(R.id.sweep_retake);
        mPlot = findViewById(R.id.sweep_plot);
        mText = findViewById(R.id.sweep_text);

        mDistSpkr = getIntent().getFloatExtra("distSpkr", 0.0f)/100;
        mHeightDiff = getIntent().getFloatExtra("heightDiff", 0.0f)/100;
        mXYDiff = Math.sqrt(Math.pow(mDistSpkr, 2) - Math.pow(mHeightDiff, 2));
        mIP = getIntent().getStringExtra("ip");

        Log.d(LOG_TAG, "Passed in a distance to speaker of " + Double.toString(mDistSpkr) +
                " and a height difference of " + Double.toString(mHeightDiff));

        mText.setText(String.format("%s%s", getResources().getString(R.string.text_mic),
                Integer.toString(mIRidx + 1)));
        mBut.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                sweepButton(v);
            }
        });
        mRetake.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                if (mIRidx > 0) {
                    mIRidx--;
                    sweepButton(v);
                }
            }
        });
        mPlot.getLegend().setVisible(false);
        //PanZoom.attach(mPlot, PanZoom.Pan.BOTH, PanZoom.Zoom.STRETCH_HORIZONTAL);
        PanZoom.attach(mPlot);

        mIRs = new double[NUM_MICS][];
        mIRpeaks = new int[NUM_MICS][];
        mIRDirects = new int[NUM_MICS];

    }

    int BufferElements2Rec = 1024; // want to play 2048 (2K) since 2 bytes we use only 1024
    int BytesPerElement = 2; // 4 bytes in IEEE754 format

    private void sweepButton(View v) {
        if (!mIsRecording && mIRidx < NUM_MICS) {
            startSweep();
        } else if (mIRidx >= NUM_MICS) {
            Bundle bundle = new Bundle();
            for (int i = 0; i < NUM_MICS; i++) {
                String label = "peaks" + Integer.toString(i+1);
                bundle.putIntArray(label, mIRpeaks[i]);
            }
            bundle.putString("ip", mIP);
            bundle.putInt("fs", 48000);
            Intent intent = new Intent(v.getContext(), RoomShapeActivity.class);
            intent.putExtras(bundle);
            startActivity(intent);
        }
    }
    private void startSweep() {

        // UI
        mBut.setText(R.string.button_sweeping);
        mBut.setEnabled(false);
        mText.setText(R.string.text_sweep);

        mRecorder = new AudioRecord(MediaRecorder.AudioSource.MIC,
                FS, RECORDER_CHANNELS,
                RECORDER_AUDIO_ENCODING, BufferElements2Rec * BytesPerElement);
        mPlayer = MediaPlayer.create(getApplicationContext(), R.raw.sweep);

        mPlayer.setOnCompletionListener(new MediaPlayer.OnCompletionListener() {
            @Override
            public void onCompletion(MediaPlayer mediaPlayer) {
                Handler handler = new Handler();
                handler.postDelayed(new Runnable() {
                    @Override
                    public void run() {
                            stopRecording();
                    }
                }, 1000);
            }
        });
        mRecorder.startRecording();
        mIsRecording = true;
        mRecordingThread = new Thread(new Runnable() {
            public void run() {
                writeAudioDataToBuffer();
            }
        }, "AudioRecorder Thread");
        mRecordingThread.start();
        mPlayer.start();
    }

    private void writeAudioDataToBuffer() {
        // Write the output audio in byte
        float sData[] = new float[BufferElements2Rec];
        mAudioBuf = FloatBuffer.allocate(FS * 15);
        while (mIsRecording) {
            mRecorder.read(sData, 0, BufferElements2Rec, AudioRecord.READ_BLOCKING);
            try {
                mAudioBuf.put(sData);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private void stopRecording() {

        // UI
        mBut.setText(R.string.button_calculating);
        mText.setText(R.string.text_ir_calc);

        // stops the recording activity
        if (null != mRecorder) {
            mIsRecording = false;
            //mRecorder.stop();
            mRecorder.release();
            mRecorder = null;
            mRecordingThread = null;
        }

        Thread correlationThread = new Thread(new Runnable() {
            @Override
            public void run() {
                doCorrelation();
                mIRidx++;
            }
        });
        correlationThread.start();
    }

    private void doCorrelation() {

        // Step 1: Convert recorded data and reference data to zero-padded double arrays
        Log.d(LOG_TAG, "Begin correlation");
        // Get raw resource from android and find length
        WavFile wavFile;
        try {
            wavFile = WavFile.openWavFile(getResources().openRawResourceFd(R.raw.sweep));
        } catch (IOException | WavFileException e) {
            e.printStackTrace();
            return;
        }
        Log.d(LOG_TAG, "Wavfile has " + Integer.toString(wavFile.getNumChannels()) +
                " channels and a total of " + Long.toString(wavFile.getNumFrames()) + " frames.");

        // Find max length and allocate arrays
        int signalLen = Math.max(mAudioBuf.position(), (int)wavFile.getNumFrames());
        double[] refData = new double[signalLen];
        double[] recData = new double[signalLen];

        // Get reference data
        double[] doubleBuffer = new double[1024];
        int doublesRead;
        int i = 0;
        try {
            while ((doublesRead = wavFile.readFrames(doubleBuffer, 1024)) > 0) {
                for (int j = 0; j < doublesRead; j++) {
                    refData[i + j] = doubleBuffer[j];
                }
                i += doublesRead;
            }
        } catch (IOException | WavFileException e) {
            e.printStackTrace();
        }
        for (; i < refData.length; i++) {
            refData[i] = 0.0f;
        }

        // Get recorded data
        for (i = 0; i < mAudioBuf.position(); i++) {
            recData[i] = mAudioBuf.get(i);
        }
        for (; i < recData.length; i++) {
            recData[i] = 0.0f;
        }
        Log.d(LOG_TAG, "All data processed. Padded reference data size " +
                Long.toString(wavFile.getNumFrames()) + " and recorded data size " +
                Integer.toString(mAudioBuf.position()) + " both to padded size " +
                Integer.toString(signalLen));

        // Step 2: FFT the recorded data and reference data
        DoubleFFT_1D fft = new DoubleFFT_1D(signalLen);
        fft.realForward(refData);
        double[] refFFT = refData.clone();
        fft.realForward(recData);
        double[] recFFT = recData;
        double[] sumFFT = new double[signalLen];

        Log.d(LOG_TAG, "Forward FFTs both complete");

        // Step 3: Mulitply reference and conj(recorded) and IFFT
        double a, b, c, d, W;
        for (i = 0; i < signalLen; i+=2) {
            a = recFFT[i];
            b = -recFFT[i+1];
            c = refFFT[i];
            d = refFFT[i+1];
            sumFFT[i] = (a*c) - (b*d);
            sumFFT[i+1] = (b*c) + (a*d);
            W = 1 / (complex_mag(sumFFT[i], sumFFT[i+1]) + 1e-15F);
            sumFFT[i] = sumFFT[i]*W;
            sumFFT[i+1] = sumFFT[i+1]*W;
        }

        // Done beforehand, because doing after realInverse makes all doubles 0? I don't know.
        Log.d(LOG_TAG, "Complex multiply complete");
        int impulseOffset = (int)Math.round(distToSpkrInM(mIRidx+1) *
                (FS/SPEED_OF_SOUND));
        Log.d(LOG_TAG,"Expected impulse offset: " + Integer.toString(impulseOffset));
        fft.realInverse(sumFFT, true);
        double[] impulseFull = new double[signalLen];
        for (i = 0; i < signalLen; i++) {
            impulseFull[i] = sumFFT[signalLen-i-1];
        }
        double maxVal = 0.0;
        mIRDirects[mIRidx] = impulseOffset;
        for (i = impulseOffset; i < impulseFull.length; i++) {
            if (Math.abs(impulseFull[i]) > Math.abs(impulseFull[mIRDirects[mIRidx]])) {
                mIRDirects[mIRidx] = i;
                maxVal = Math.abs(impulseFull[i]);
            }
        }
        if (mIRDirects[mIRidx] > 20000) {
            Log.w(LOG_TAG, "Unexpected direct signal way out of normal range! Sample index: " +
            Integer.toString(mIRDirects[mIRidx]));
        }
        int startSamplesToSkip = mIRDirects[mIRidx] - impulseOffset;
        double scale = 1.0/maxVal;
        mIRDirects[mIRidx] = impulseOffset;
        mIRs[mIRidx] = new double[signalLen-startSamplesToSkip];
        for (i = 0; i < mIRs[mIRidx].length; i++) {
            mIRs[mIRidx][i] = Math.exp(Math.abs(impulseFull[i+startSamplesToSkip])*scale);
        }
        Log.d(LOG_TAG, "Inverse FFT complete. Finding peaks.");

        // Step 4. Find peaks
        int window_half = Math.round(WINDOW_SIZE/2);
        // peaks, in samples from the highest peak in the signal
        ArrayList<Integer> peaks = new ArrayList<>();
        peaks.add(mIRDirects[mIRidx]);
        int maxSample = mIRDirects[mIRidx];
//        for (i = maxSample/window_half; i < (maxSample+WINDOW_POST)/window_half; i++) {
//            int window_start = i*window_half;
        for (int window_start = maxSample;
             window_start < maxSample+WINDOW_POST;
             window_start += window_half) {
            double[] window = new double[WINDOW_SIZE];
            System.arraycopy(mIRs[mIRidx], window_start, window, 0, WINDOW_SIZE);
            int maxIdx = 0;
            maxVal = NEGATIVE_INFINITY;
            for (int j = 0; j < window.length; j++) {
                if (window[j] > maxVal) {
                    maxIdx = j;
                    maxVal = window[j];
                }
            }
            if (maxIdx > (WINDOW_SIZE/4) && maxIdx <= 3*(WINDOW_SIZE/4)
                    && maxVal > PEAK_THRESHOLD) {
                peaks.add(window_start + maxIdx);
            }
        }
        mIRpeaks[mIRidx] = new int[peaks.size()];
        for (i = 0; i < mIRpeaks[mIRidx].length; i++) {
            mIRpeaks[mIRidx][i] = peaks.get(i);
        }
        Log.d(LOG_TAG, "Peaks found, graphing...");

        // Step 5. Graph it
        Integer[] IRXvals = new Integer[IR_LEN];
        Double[] IRYvals = new Double[IR_LEN];
        for (i = 0; i < IRYvals.length; i++) {
            IRXvals[i] = i*IR_HSCALE;
            IRYvals[i] = mIRs[mIRidx][i];
        }

        XYSeries IRseries = new SimpleXYSeries(Arrays.asList(IRYvals),
                SimpleXYSeries.ArrayFormat.Y_VALS_ONLY, "Impulse Response");
        LineAndPointFormatter IRformatter = new LineAndPointFormatter(getApplicationContext(),
                R.xml.line_format);

        Integer[] peakXvals = new Integer[peaks.size()];
        Double[] peakYvals = new Double[peaks.size()];
        for (i = 0; i < peaks.size(); i++) {
            peakXvals[i] = peaks.get(i);
            peakYvals[i] = mIRs[mIRidx][peaks.get(i)];
        }
        XYSeries peaksSeries = new SimpleXYSeries(Arrays.asList(peakXvals),
                Arrays.asList(peakYvals), "Peaks");
        LineAndPointFormatter peaksFormatter = new LineAndPointFormatter(getApplicationContext(),
                R.xml.point_format);
        mPlot.clear();
        mPlot.addSeries(IRseries, IRformatter);
        mPlot.addSeries(peaksSeries, peaksFormatter);
        mPlot.redraw();
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if (mIRidx <= NUM_MICS-1) {
                    mText.setText(String.format("%s%s", getResources().getString(R.string.text_mic),
                            Integer.toString(mIRidx + 1)));
                    mBut.setText(R.string.button_go);
                } else {
                    mText.setText(R.string.text_calculate_room_shape);
                    mBut.setText(R.string.button_graph);
                }
                mBut.setEnabled(true);
                mRetake.setEnabled(true);
            }
        });
        Log.d(LOG_TAG, "Graph shown");
    }

    float[] fft_abs(float[] fft_in) {
        float[] abs = new float[fft_in.length/2];
        for (int i = 0; i < abs.length; i++) {
            abs[i] = complex_mag(fft_in[2*i], fft_in[(2*i)+1]);
        }
        return abs;
    }

    double[] fft_abs(double[] fft_in) {
        double[] abs = new double[fft_in.length/2];
        for (int i = 0; i < abs.length; i++) {
            abs[i] = complex_mag(fft_in[2*i], fft_in[(2*i)+1]);
        }
        return abs;
    }

    float complex_mag(float real, float imag) {
        return (float)Math.sqrt(Math.pow(real, 2) +
                                Math.pow(imag, 2));
    }

    double complex_mag(double real, double imag) {
        return Math.sqrt(Math.pow(real, 2) +
                         Math.pow(imag, 2));
    }

    double distToSpkrInM(int micIdx) {
        /*
        Spiral: if you're facing mic from speaker
        Left 10, away 20, right 30, towards 45
         */

        double l;
        double w;
        double xy;
        switch (micIdx) {
            case 1:
                return mDistSpkr;
            case 2:
                l = mDistSpkr;
                w = 0.1;
                xy = Math.sqrt(Math.pow(l, 2) + Math.pow(w, 2));
                return Math.sqrt(Math.pow(xy, 2) + Math.pow(mHeightDiff, 2));
                //return Math.sqrt(Math.pow(xy, 2) + Math.pow(Math.abs(mHeightDiff-57), 2));
            case 3:
                l = mXYDiff + 0.2;
                w = 0.1;
                xy = Math.sqrt(Math.pow(l, 2) + Math.pow(w, 2));
                return Math.sqrt(Math.pow(xy, 2) + Math.pow(mHeightDiff, 2));
            case 4:
                l = mXYDiff + 0.2;
                w = 0.2;
                xy = Math.sqrt(Math.pow(l, 2) + Math.pow(w, 2));
                return Math.sqrt(Math.pow(xy, 2) + Math.pow(mHeightDiff, 2));
                //return Math.sqrt(Math.pow(xy, 2) + Math.pow(Math.abs(mHeightDiff-57), 2));
            case 5:
                l = mXYDiff - 0.25;
                w = 0.2;
                xy = Math.sqrt(Math.pow(l, 2) + Math.pow(w, 2));
                return Math.sqrt(Math.pow(xy, 2) + Math.pow(mHeightDiff, 2));
            default:
                Log.w(LOG_TAG, "distToSpkrInM called with invalid mic index!");
                return 0.0;
        }
    }
}
