/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.help;

import com.google.gson.Gson;
import org.apache.commons.io.IOUtils;
import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.DefaultHttpClient;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.List;

public class ForumAPIUtils {
    /**
     * How we post to the forum
     */
    final private static String ACCESS_TOKEN = "access_token=";

    public static List<String> getPostedTools(String forumKey) {
        Gson gson = new Gson();
        List<String> output = new ArrayList<String>();

        String text = httpGet(HelpConstants.GATK_FORUM_API_URL + "categories.json?CategoryIdentifier=tool-bulletin&page=1-100000&" + ACCESS_TOKEN + forumKey);

        APIQuery details = gson.fromJson(text, APIQuery.class);
        ForumDiscussion[] discussions = details.Discussions;

        for (ForumDiscussion post : discussions) {
            output.add(post.Name);
        }

        /*
        System.out.println(details.isJsonArray());
        System.out.println(details.isJsonNull());
        System.out.println(details.isJsonObject());
        System.out.println(details.isJsonPrimitive());

        JsonArray posted = details.getAsJsonPrimitive().get("Discussions").getAsJsonArray();

        for ( JsonElement post : posted ) {
            output.add( post.getAsJsonObject().get("Name").getAsString());
        }
        */
        return output;
    }


    private static String httpGet(String urlStr) {
        String output = "";

        try {

            DefaultHttpClient httpClient = new DefaultHttpClient();
            HttpGet getRequest = new HttpGet(urlStr);
            getRequest.addHeader("accept", "application/json");

            HttpResponse response = httpClient.execute(getRequest);

            if (response.getStatusLine().getStatusCode() != 200) {
                throw new RuntimeException("Failed : HTTP error code : "
                        + response.getStatusLine().getStatusCode());
            }

            output = IOUtils.toString(response.getEntity().getContent());

            httpClient.getConnectionManager().shutdown();

        } catch (ClientProtocolException e) {

            e.printStackTrace();

        } catch (IOException e) {

            e.printStackTrace();
        }
        return output;
    }

    private static String httpPost(String data, String URL) {
        try {

            DefaultHttpClient httpClient = new DefaultHttpClient();
            HttpPost postRequest = new HttpPost(URL);

            StringEntity input = new StringEntity(data);
            input.setContentType("application/json");
            postRequest.setEntity(input);

            HttpResponse response = httpClient.execute(postRequest);

            if (response.getStatusLine().getStatusCode() != 200) {
                throw new RuntimeException("Failed : HTTP error code : "
                        + response.getStatusLine().getStatusCode());
            }

            BufferedReader br = new BufferedReader(
                    new InputStreamReader((response.getEntity().getContent())));

            String output = "";
            String line;
            System.out.println("Output from Server .... \n");
            while ((line = br.readLine()) != null) {
                output += (line + '\n');
                System.out.println(line);
            }

            br.close();
            httpClient.getConnectionManager().shutdown();
            return output;

        } catch (MalformedURLException e) {

            e.printStackTrace();

        } catch (IOException e) {

            e.printStackTrace();

        }
        return null;
    }

    public static void postToForum(GATKDocWorkUnit tool, final String forumKey) {


        ForumDiscussion post = new ForumDiscussion(tool);

        Gson gson = new Gson();

        String data = gson.toJson(post.getPostData());
        httpPost(data, HelpConstants.GATK_FORUM_API_URL + "post/discussion.json?" + ACCESS_TOKEN + forumKey);


    }

    class APIQuery {
        ForumDiscussion[] Discussions;

        public APIQuery() {}
    }

}
