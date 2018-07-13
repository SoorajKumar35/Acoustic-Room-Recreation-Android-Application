package com.bitbucket.arrc.arrc_acousticroomrecreation;

import android.app.Activity;
import android.content.Intent;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;

import org.json.JSONException;
import org.json.JSONObject;

import java.io.DataOutputStream;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;

public class RoomShapeActivity extends Activity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_room_shape);
        final TextView text = findViewById(R.id.roomshape_text);
        Button but = findViewById(R.id.roomshape_restart);

        but.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                startActivity(new Intent(v.getContext(), MainActivity.class));
            }
        });

        final int[] peaks1 = getIntent().getIntArrayExtra("peaks1");
        final int[] peaks2 = getIntent().getIntArrayExtra("peaks2");
        final int[] peaks3 = getIntent().getIntArrayExtra("peaks3");
        final int[] peaks4 = getIntent().getIntArrayExtra("peaks4");
        final int[] peaks5 = getIntent().getIntArrayExtra("peaks5");
        final String ip = getIntent().getStringExtra("ip");
        final int fs = getIntent().getIntExtra("fs", 48000);
        Thread thread = new Thread(new Runnable() {
            @Override
            public void run() {
                final String res = roomshape(fs, peaks1, peaks2, peaks3, peaks4, peaks5);
                try {
                    URL url = new URL(ip);
                    HttpURLConnection conn = (HttpURLConnection) url.openConnection();
                    conn.setRequestMethod("POST");
                    conn.setRequestProperty("Content-Type", "application/json;charset=UTF-8");
                    conn.setRequestProperty("Accept", "application/json");
                    conn.setDoOutput(true);
                    conn.setDoInput(true);

                    DataOutputStream os = new DataOutputStream(conn.getOutputStream());
                    os.writeBytes(res);

                    os.flush();
                    os.close();

                    Log.i("STATUS", String.valueOf(conn.getResponseCode()));
                    Log.i("MSG", conn.getResponseMessage());

                    conn.disconnect();

                } catch (IOException e) {
                    e.printStackTrace();
                }
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        text.setText(res);
                        Log.d("ARRC", res);
                    }
                });
            }
        });
        thread.start();
    }

    public native String roomshape(int fs, int[] peaks1, int[] peaks2,
                                   int[] peaks3, int[] peaks4, int[] peaks5);

    static {
        System.loadLibrary("roomshape");
    }
}
